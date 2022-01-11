close all
clearvars

%% Virtual patients

%loading patient parameters
load VirtualPatients_11_01_22
VirtualP=[VP];

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
    
    p.dosePAC1=75;
    p.doseTRAIL=3*VirtualP(j,5);
    
    %  Update model parameters
    [p.(ParameterNames{1}),p.(ParameterNames{2}),p.(ParameterNames{3})] = C{:}; %update the parameters for this run 

     % determine optimal dose for patient 
    FObjective=@(x)ProtocolOptimizerObjective(x,p);
    opts = optimoptions('ga','EliteCount',15,'MaxStallGenerations',4,'FunctionTolerance',1e-6, 'MaxGenerations',50,'Display','iter','PopulationSize',60);
    [x,fval] = ga(FObjective, nvars,[],[],[],[],LB,UB,[],IntCon,opts);   
    
    %Store the optimal dosing regime 
    OptimalDose(j,:) = x;
    Results{j}.OptimalDose=OptimalDose;
    Results{j}.fval=fval;
    
    %simulate model on optimal dosage and store
    [sol] = modelsimulator_multidose(x,p)   
    SolMat{j} = sol;
    save('ResultsDec17','Results','SolMat')
     
    
end

%% figure plotting heatmap of optimal dosages

%load('ResultsDec17')
%OptimalDose = UpdatedResults;

sumX = [247,252,240
224,243,219
204,235,197
168,221,181
123,204,196
78,179,211
43,140,190
8,104,172
8,64,129]/255;

%Heatmap of patients optimal dosage for PAC-1
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

% heatmap of patients optimal dose for TAIL
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

%% Plotting virtual patients but segregated by High PAC-1 elimination and Normal PAC-1 elimination 

% loading IDs for patients in High PAC-1 elimianation cohort (cohort 1) and
% normal PAC-1 elimination (cohort 2)
load('cohort1_id.mat')
cohort2_id = [];
for ii = 1:400
    if isempty(find(cohort1_id == ii))==1
        cohort2_id = [cohort2_id;ii];  
    end
end

sumX = [247,252,240
224,243,219
204,235,197
168,221,181
123,204,196
78,179,211
43,140,190
8,104,172
8,64,129]/255;

% plotting heatmap of optimal PAC-1 dosage where patients are plotted in
% their subcohort
figure
h = heatmap(OptimalDose([cohort2_id;cohort1_id],1:11),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')
ylabel('Virtual Patients')
xlabel('Time (weeks)')
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
title('PAC-1')

% plotting heatmap of optimal TRAIL dosage where patients are plotted in
% their subcohort
figure
h = heatmap(OptimalDose([cohort2_id;cohort1_id],12:22),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')
ylabel('Virtual Patients')
xlabel('Time (weeks)')
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
title('TRAIL')


%% plotting patient dynamics under optimal protocol 

% reloading patient variable dynamics
for i = 1:400
    
    sol = SolMat{i};
    time = linspace(0,p.tf,1000);
    solind_C(i,:) = deval(sol,time,1);
    solind_PAC1(i,:) = deval(sol,time,3);
    solind_TRAIL(i,:) = deval(sol,time,4);
    
    CumulativeTumourBurden_mat_cohort2(i,:) = trapz(time,solind_C(i,:));
    
end

%plotting the log cumulative tumour burden
figure
histogram(log(real(CumulativeTumourBurden_mat_cohort2)),200,'FaceColor',[10 147 150]/255,'EdgeColor','none')
xlabel('log(Cumulative Tumour Burden)')
ylabel('Frequency') 
set(gca,'FontSize',18)
box off

% plotting the number of tumour cells for each patient
fig1 = figure
hold on
plot(time,solind_C,'Color',[0.5 0.5 0.5])
xlabel('Time')
ylabel('Cancer cells C, log(cells)')
set(gca,'FontSize',18)
set(gca,'yscale','log')

% plotting the mean and std for the cohort
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

% plotting mean and std for cohort PAC-1 
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
ylabel('PAC-1 plasma ({\mu}M)')
set(gca,'FontSize',18)
ylim([0 2500])

% plotting mean and std for cohort TRAIL
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
ylabel('TRAIL plasma (ng/ml)')
set(gca,'FontSize',18)
ylim([0 400])


%% calculating day of recovery and percentage recovered
recovery = [];

%calculating if patient pcount recovered, i.e. less than 1 cell remaining
for pcount = 1:400
    if isempty(find(solind_C(pcount,:)<1))==0
        loc_recovery = find(solind_C(pcount,:)<1,1);
        recovery = [recovery;time(loc_recovery)];
    end
end

%plotting a histogram for patient recovery times
figure
h = histogram(recovery)
ylabel('Frequency')
xlabel('Day of recovery')
grid on 
box off
set(gca,'FontSize',18)

%plotting the cumulative percentage of virtual patients that recover
figure
plot(h.BinEdges(1:10)+0.5,cumsum(h.Values)/400*100,'o:','Color',[0.51 0.87 0.85],'LineWidth',2)
box off
xlabel('Time (days)')
ylabel('Percentage recovered')
set(gca,'FontSize',18)

