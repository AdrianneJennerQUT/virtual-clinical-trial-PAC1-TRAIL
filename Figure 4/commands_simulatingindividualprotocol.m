% simulating patients under individualised protocol

close all
clearvars

%% Virtual patients

%loading patient parameters
load VP.mat
load VP_May31.mat
VirtualP=[VP; VP_May31];

%% Drug PK parameters

%general data
p.C0 = 5000;   %initial number of cells
p.D0=0;     %initial number of dead cells
tspan=[0 12]; %0 to 12 hours 

%KGN Gompertz model
p.r = 0.052696;
p.K = 1745423909.2602;

%mean bodyweight
p.BW=70;%kg

%PAC_1 with TRAIL data: 
IC50_PAC1_TRAIL=40.3135;
Imax_PAC1_TRAIL=1.0304;
h_PAC1_TRAIL=0.6829;
%PAC-1 data
delta_PAC1=0.02208;

%Trail data
IC50_TRAIL=15.1509;
Imax_TRAIL=0.5278;
h_TRAIL=1.0846;
delta_TRAIL=0.027865;

%drug parameters (drug 1 is PAC-1 and drug 2 is TRAIL)
p.IC50_1=IC50_PAC1_TRAIL;
p.gamma_1=h_PAC1_TRAIL;
p.delta_1=delta_PAC1;
p.Imax_1=Imax_PAC1_TRAIL;

p.IC50_2=IC50_TRAIL;
p.gamma_2=h_TRAIL;
p.delta_2=delta_TRAIL;
p.Imax_2=Imax_TRAIL;

%PAC-1 PKs
p.Cl1=24*574.7/1000; %L/kg/day
p.Vd1=0.7; %(L/kg) vd of central compartment
p.ke1=p.Cl1/p.Vd1;

%TRAIL PKs
p.Cl2=116; %L/day
p.Vd2=4.28; %L
p.ke2=p.Cl2/p.Vd2;

%Check this!!!
p.psi=1; %factor for synergy/antagonism/addition if there are two drugs (psi=1 no effect)
p.xi=p.IC50_1/p.IC50_2;

% Time interval
p.t_treat=2.5*30;
p.t_after=0.5*1;
p.tf = p.t_treat+p.t_after; %Simulation over 2.5 months of treatment + 0.5 after treatment
p.totaltime = [0 p.tf];

%PAC-1 dosing 
p.StartTimePAC1 = 0;
p.OffsetPAC1 = 7; %Can be dosed every week

%TRAIL dosing
p.StartTimeTrail = 0;
p.OffsetTrail = 7; %Can be dosed every week

%Total number of administrations
p.AdministrationTimes=0:p.OffsetPAC1:p.t_treat;
p.AdminNumberPAC1 = numel(p.AdministrationTimes);%Number of PAC-1 doses over 2.5 months
p.AdminNumberTrail = numel(p.AdministrationTimes);%Number of TRAIL  doses over 2.5 months

%Input for optimizer and storage vectors
nvars = p.AdminNumberPAC1 + p.AdminNumberTrail;
LB = 0.*ones(1,nvars);   
UB = 10.*ones(1,nvars); 
IntCon = 1:nvars; %The condition that ensures that the optimial solution enforces integer multiple doses of the baseline dose

load('UpdatedResults')
OptimalDose = UpdatedResults; %Matrix to store the optimal doses for each patient


for j=1:size(VirtualP,1) %For each patient, find the optimal dosing regime
    
    disp(j)
    
    VPload = VirtualP(j,[1,2,3,4]); %load the varied parameters from the virtual patient and update VPx.
    %calculate ke for each drug from Cl and Vd
    p.ke1=VPload(2)/VPload(1);
    p.ke2=VPload(4)/VPload(3);
    %need the ke for each drug and Vd for PAC1 (to calculate dose)
    Vd1=VPload(1);
    VParam=[p.ke1 p.Vd1 p.ke2];
    C = {VParam(1), ... 
    VParam(2), ...
    VParam(3)};
    ParameterNames = {'ke1','Vd1','ke2'}; %The name of parameters to be varied
    
    p.dosePAC1=75;
    p.doseTRAIL=3*VirtualP(j,5);
    
    %  Update model parameters
    [p.(ParameterNames{1}),p.(ParameterNames{2}),p.(ParameterNames{3})] = C{:}; %update the parameters for this run 

    %Store the optimal dosing regime 
    x=OptimalDose(j,:);
    
    %simulate model on optimal dosage and store
    [sol] = modelsimulator_multidose(x,p);   
    SolMat{j} = sol;
    save('ResultsDec17','SolMat')
     
    
end


