close all
clearvars

%% Virtual patients

% %uncomment to regenerate a new cohort
% num_patients=100;
% VP=GeneratingPatients(num_patients);

load VP.mat
load VP_May31.mat
VirtualP=[VP; VP_May31];

%% Drug PK parameters

%general data
p.C0 = 5000;   %initial number of cells
p.D0=0;     %initial number of dead cells
tspan=[0 12]; %0 to 12 hours 
%p.r = 0.27073;%Gompertz growth %0.3493; %proliferation rate obtained in vitro (per day)
%p.K = 1607.5956*1E6;
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
OptimalDose = zeros(size(VirtualP,1),nvars); %Matrix to store the optimal doses for each patient

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
    
    p.dosePAC1=3*VirtualP(j,5);
    p.doseTRAIL=75;
    
    %  Update model parameters
    [p.(ParameterNames{1}),p.(ParameterNames{2}),p.(ParameterNames{3})] = C{:}; %update the parameters for this run 

    tic 
    FObjective=@(x)ProtocolOptimizerObjective(x,p);
    opts = optimoptions('ga','EliteCount',15,'MaxStallGenerations',4,'FunctionTolerance',1e-6, 'MaxGenerations',50,'Display','iter','PopulationSize',60);
    [x,fval] = ga(FObjective, nvars,[],[],[],[],LB,UB,[],IntCon,opts);   
    
    %Store the optimal dosing regime ****
    OptimalDose(j,:) = x;
    Results{j}.OptimalDose=OptimalDose;
    Results{j}.fval=fval;
    
    %simulate model on optimal dosage and store
    [sol] = modelsimulator_multidose(x,p)   
    SolMat{j} = sol;
    save('Results_Nov8','Results','SolMat')
    
    
 
    toc
    
end

%% figure plotting heatmap of optimal dosages

sumX = [0 18 25
    0 95 115
    10 147 150
    148 210 189
    233 216 166
    238 155 0
    202 103 2
    187 62 3
    174 32 18
    155 34 38]/255;

sumf = [255,247,243
253,224,221
252,197,192
250,159,181
247,104,161
221,52,151
174,1,126
122,1,119
73,0,106]/255;

figure
h = heatmap(OptimalDose(:,1:11),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')
ylabel('Virtual Patients')
xlabel('Weekly dose')
title('PAC-1')
set(gca,'FontSize',16)
h.XDisplayLabels(1) = {"0"}
h.XDisplayLabels(2) = {"1"}
h.XDisplayLabels(3) = {"2"}
h.XDisplayLabels(4) = {"3"}
h.XDisplayLabels(5) = {"4"}
h.XDisplayLabels(6) = {"5"}
h.XDisplayLabels(7) = {"6"}
h.XDisplayLabels(8) = {"7"}
h.XDisplayLabels(9) = {"8"}
h.XDisplayLabels(10) = {"9"}
h.XDisplayLabels(11) = {"10"}

h.YDisplayLabels(1:400) = {""}

s = struct(h);
s.XAxis.TickLabelRotation = 0;  % vertical


figure
h = heatmap(OptimalDose(:,12:22),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')
ylabel('Virtual Patients')
xlabel('Weekly dose')
title('TRAIL')
set(gca,'FontSize',18)
h.XDisplayLabels(1) = {"0"}
h.XDisplayLabels(2) = {"1"}
h.XDisplayLabels(3) = {"2"}
h.XDisplayLabels(4) = {"3"}
h.XDisplayLabels(5) = {"4"}
h.XDisplayLabels(6) = {"5"}
h.XDisplayLabels(7) = {"6"}
h.XDisplayLabels(8) = {"7"}
h.XDisplayLabels(9) = {"8"}
h.XDisplayLabels(10) = {"9"}
h.XDisplayLabels(11) = {"10"}

h.YDisplayLabels(1:400) = {""}

s = struct(h);
s.XAxis.TickLabelRotation = 0;  % vertical

%% plotting patient dynamics under optimal protocol 

for i = 1:400
    
    sol = SolMat{i};
    time = linspace(0,p.tf,1000);
    solind_C(i,:) = deval(sol,time,1);
    solind_D(i,:) = deval(sol,time,2);
    solind_PAC1(i,:) = deval(sol,time,3);
    solind_TRAIL(i,:) = deval(sol,time,4);
    
    CumulativeTumourBurden_mat_cohort2(i,:) = trapz(time,solind_C(i,:));
    
end


figure
histogram(log(real(CumulativeTumourBurden_mat_cohort2)),200,'FaceColor',[10 147 150]/255,'EdgeColor','none')
xlabel('log(Cumulative Tumour Burden)')
ylabel('Frequency') 
set(gca,'FontSize',18)
box off


fig1 = figure
hold on
plot(time,solind_C,'Color',[0.5 0.5 0.5])
xlabel('Time')
ylabel('Cancer cells C, log(cells)')
set(gca,'FontSize',18)
set(gca,'yscale','log')

fig1 = figure
hold on
options.handle = fig1;
options.color_area = [0 95 115]/255;
options.color_line =  [0 95 115]/255;
options.alpha = 0.3;
options.error = 'std';
options.line_width = 2;
options.x_axis = time;
plot_areaerrorbar(real(solind_C),options)
xlabel('Time')
ylabel('Cancer cells C, cells')
set(gca,'FontSize',18)
ylim([0 2e6])

fig2 = figure
hold on
options.handle = fig2;
options.color_area = [252,141,98]/255;
options.color_line = [252,141,98]/255;
options.alpha = 0.3;
options.error = 'std';
options.line_width = 2;
options.x_axis = time;
plot_areaerrorbar(real(solind_D),options)
xlabel('Time')
ylabel('Dead cells D, cells')
set(gca,'FontSize',18)

fig1 = figure
hold on
options.handle = fig1;
options.color_area = [202 103 2]/255;
options.color_line = [202 103 2]/255;
options.alpha = 0.3;
options.error = 'std';
options.line_width = 2;
options.x_axis = time;
plot_areaerrorbar(solind_PAC1,options)
xlabel('Time')
ylabel('PAC-1')
set(gca,'FontSize',18)
ylim([0 2500])


fig1 = figure
hold on
options.handle = fig1;
options.color_area = [155 34 38]/255;
options.color_line = [155 34 38]/255;
options.alpha = 0.3;
options.error = 'std';
options.line_width = 2;
options.x_axis = time;
plot_areaerrorbar(solind_TRAIL,options)
xlabel('Time')
ylabel('TRAIL')
set(gca,'FontSize',18)
ylim([0 400])




