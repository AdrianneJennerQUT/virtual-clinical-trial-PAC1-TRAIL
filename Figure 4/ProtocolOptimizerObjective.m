%% Fitness function
function [Obj] = ProtocolOptimizerObjective(x,PA)   

%PAC-1
%Load the dosing strategy for PAC-1
PAC1Dose = x(1:PA.AdminNumberPAC1);
PA.DosePAC1 = PA.dosePAC1.*(PAC1Dose); %Amount of IL2 administered (input from the optimizer)

%TRAIL 
%Load the dosing strategy for TRAIL
TrailDose = x(PA.AdminNumberPAC1+1:PA.AdminNumberTrail+PA.AdminNumberPAC1); 
PA.DoseTrail = PA.doseTRAIL.*(TrailDose); %Amount of IL7 administered (input from the optimizer)

%Calculate the total dose given:
TotalDose = dot(PA.DoseTrail,ones(1,PA.AdminNumberTrail)) + dot(PA.DosePAC1,ones(1,PA.AdminNumberPAC1)); 

%Initial conditions
PA.IC = [PA.C0;PA.D0;PA.DosePAC1(1);PA.DoseTrail(1)];
        
% Solve system
[solTreat] = SimulateComboDrugModel(PA);

%Objective function: Area under the tumour curve+area under treatment curve
TimeSeries = linspace(0,PA.tf,1001); %Create 1000 evenly space points inside the time domain
EvalSol = deval(solTreat,TimeSeries,1); %Evaluate the tumour solution at the collocation points
CumulativeTumourBurden = trapz(TimeSeries,EvalSol); %Calculate the tumour AUC

%The objective function we are trying to minimize.
Obj = CumulativeTumourBurden + TotalDose;

%% Distributed Immunity  Solver

function [sol] = SimulateComboDrugModel(PA)
    
options=[];

if numel(PA.AdministrationTimes) == 1
    sol=ode15s(@(t,y)odefun(t,y,PA),[PA.AdministrationTimes(1) PA.tf],PA.IC,options);
else
    %Call the solver to calculate solution over initial administration
    sol=ode15s(@(t,y)odefun(t,y,PA),[0 PA.AdministrationTimes(2)],PA.IC,options);
    if numel(PA.AdministrationTimes) > 2
        %Extend the solution for each subsequent administration
        for nn=3:1:numel(PA.AdministrationTimes)
            %Set ICs to be final solution in last interval
            PA.IC=sol.y(:,end);
            PA.IC(1)=PA.IC(1);
            PA.IC(2)=PA.IC(2);
            PA.IC(3)=PA.IC(3)+PA.DosePAC1(nn-1);
            PA.IC(4)=PA.IC(4)+PA.DoseTrail(nn-1);
            sol=odextend(sol,@(t,y)odefun(t,y,PA),PA.AdministrationTimes(nn),PA.IC,options);
        end
    %       Calculate solution from end of last admin to final time
        PA.IC=sol.y(:,end);
        sol=odextend(sol,@(t,y)odefun(t,y,PA),PA.tf,PA.IC,options);
    else
        %Set ICs to be final solution in last interval
        PA.IC=sol.y(:,end);
        %Calculate solution from end of last admin to final time
        sol=odextend(sol,@(t,y)odefun(t,y,PA),PA.tf,PA.IC,options);
    end
end 

function dydt=odefun(t,y,PA)
    
C = y(1);
D = y(2);
Cp1 = y(3);
Cp2 = y(4);    

exp1 = Cp1.^PA.gamma_1./(PA.psi*PA.IC50_1).^PA.gamma_1;
exp2 = (PA.xi*Cp2).^PA.gamma_2./(PA.psi*PA.IC50_1).^PA.gamma_2;

E = (PA.Imax_1.*exp1+PA.Imax_2.*exp2+(PA.Imax_1+PA.Imax_2-PA.Imax_1*PA.Imax_2)*exp1.*exp2)./(exp1+exp2+exp1.*exp2+1);%E fct for two drug 

dC = (C.*PA.r.*log(PA.K./PA.C0).*exp(-1.*PA.r.*t))-C*E;%PA.r*C-C*E;
dD = C*E;
dCp1 = -PA.ke1*Cp1;
dCp2 = -PA.ke2*Cp2; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble vector of RHS to solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dydt=[dC
      dD
      dCp1
      dCp2];
end

end

end