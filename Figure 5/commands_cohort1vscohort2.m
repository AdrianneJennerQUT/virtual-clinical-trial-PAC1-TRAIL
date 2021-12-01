% In this code, the dynamics of the cohorts for High PAC-1 elimination and Normal PAC-1
% elimination virtual patients are examined in more detail

%load virtual patient parameters
load VP.mat
load VP_May31.mat
VirtualP=[VP; VP_May31];

%load subcohort IDs
load('cohort1_id.mat')
cohort2_id = [];
for ii = 1:400
    if isempty(find(cohort1_id == ii))==1
        cohort2_id = [cohort2_id;ii];  
    end
end

cohort1_parameters = VirtualP(cohort1_id,:);
cohort2_parameters = VirtualP(cohort2_id,:);

%% plotting cohort parameters

% violin plot for Vd_PAC1 for patients in the High PAC-1 elimination cohort
% (cohort 1) and the Normal PAC-1 elimination cohort (cohort 2)
figure 
v = violinplot([[cohort1_parameters(:,1);NaN(330,1)],cohort2_parameters(:,1)]);   
ylabel('VD_{PAC1}')
title('PAC-1 volume distribution')
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;
set(gca,'xtick',[1 2],'xticklabels',{'Cohort 1','Cohort 2'})

% violin plot for CL_PAC1 for patients in the High PAC-1 elimination cohort
% (cohort 1) and the Normal PAC-1 elimination cohort (cohort 2)
figure 
v=violinplot([[cohort1_parameters(:,2);NaN(330,1)],cohort2_parameters(:,2)]);    
ylabel('CL_{PAC1}')
title('PAC-1 clearance')
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;
set(gca,'xtick',[1 2],'xticklabels',{'Cohort 1','Cohort 2'})

% violin plot for ke_PAC1 for patients in the High PAC-1 elimination cohort
% (cohort 1) and the Normal PAC-1 elimination cohort (cohort 2)
figure 
v=violinplot([[cohort1_parameters(:,2)./cohort1_parameters(:,1);NaN(330,1)],cohort2_parameters(:,2)./cohort2_parameters(:,1)]);   
ylabel('k_{ePAC1} (1/day}')
title('PAC-1 clearance rate')
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;
set(gca,'xtick',[1 2],'xticklabels',{'Cohort 1','Cohort 2'})

% violin plot for Vd_TRAIL for patients in the High PAC-1 elimination cohort
% (cohort 1) and the Normal PAC-1 elimination cohort (cohort 2)
figure 
v=violinplot([[cohort1_parameters(:,3);NaN(330,1)],cohort2_parameters(:,3)]);     
ylabel('VD_{TRAIL}') 
title('TRAIL volume distribution')
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;
set(gca,'xtick',[1 2],'xticklabels',{'Cohort 1','Cohort 2'})


% violin plot for CL_TRAIL for patients in the High PAC-1 elimination cohort
% (cohort 1) and the Normal PAC-1 elimination cohort (cohort 2)
figure 
v=violinplot([[cohort1_parameters(:,4);NaN(330,1)],cohort2_parameters(:,4)]);   
ylabel('CL_{TRAIL}')
title('TRAIL clearance')
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;
set(gca,'xtick',[1 2],'xticklabels',{'Cohort 1','Cohort 2'})

% violin plot for ke_TRAIL for patients in the High PAC-1 elimination cohort
% (cohort 1) and the Normal PAC-1 elimination cohort (cohort 2)
figure 
v=violinplot([[cohort1_parameters(:,4)./cohort1_parameters(:,3);NaN(330,1)],cohort2_parameters(:,4)./cohort2_parameters(:,3)]);   
ylabel('k_{eTRAIL} (1/day}') 
title('TRAIL clearance rate')
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;
set(gca,'xtick',[1 2],'xticklabels',{'Cohort 1','Cohort 2'})

% violin plot for BW for patients in the High PAC-1 elimination cohort
% (cohort 1) and the Normal PAC-1 elimination cohort (cohort 2)
figure 
v=violinplot([[cohort1_parameters(:,5);NaN(330,1)],cohort2_parameters(:,5)]);   
ylabel('BW (kg)')
title('Body weight')
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;
set(gca,'xtick',[1 2],'xticklabels',{'Cohort 1','Cohort 2'})


%% load results and plot cohort dynamics
load('Results_Nov8')

% Time interval
p.t_treat=2.5*30;
p.t_after=0.5*1;
p.tf = p.t_treat+p.t_after; %Simulation over 2.5 months of treatment + 0.5 after treatment
p.totaltime = [0 p.tf];

for i = 1:400
    
    sol = SolMat{i};
    time = linspace(0,p.tf,1000);
    solind_C(i,:) = deval(sol,time,1); % record total number of cancer cells
    solind_PAC1(i,:) = deval(sol,time,3); % record total PAC-1 concentration
    solind_TRAIL(i,:) = deval(sol,time,4); % record total TRAIL concentration
    
    CumulativeTumourBurden_mat(i,:) = trapz(time,solind_C(i,:)); % record cumulative tumour burden
    
end

% plot histograms comparing cumulative tumour burden for each sub cohort
figure
histogram(log(real(CumulativeTumourBurden_mat(cohort1_id,:))),10,'FaceColor',[209,105,61]/255,'EdgeColor','none')
hold on
histogram(log(real(CumulativeTumourBurden_mat(cohort2_id,:))),100,'FaceColor',[59 164 209]/255,'EdgeColor','none')
xlabel('log(Cumulative Tumour Burden)')
ylabel('Frequency') 
set(gca,'FontSize',18)
box off
legend('Cohort 1','Cohort 2')


% plot comparing number of cancer cells for each individual in each sub cohort
fig1 = figure
hold on
plot(time,solind_C(cohort1_id,:),'Color',[209,105,61]/255)
plot(time,solind_C(cohort2_id,:),'Color',[59 164 209]/255)
l1 = plot(time,solind_C(cohort1_id(1),:),'Color',[209,105,61]/255)
l2 = plot(time,solind_C(cohort2_id(1),:),'Color',[59 164 209]/255)
xlabel('Time')
ylabel('Cancer cells C, log(cells)')
set(gca,'FontSize',18)
set(gca,'yscale','log')
legend([l1 l2],{'Cohort 1','Cohort 2'})

% plot mean and std for number of cancer cells for each sub cohort
fig1 = figure
hold on
options.handle = fig1;
options.color_area = [209,105,61]/255;
options.color_line = [209,105,61]/255;
options.alpha = 0.3;
options.error = 'std';
options.line_width = 2;
options.x_axis = time;
plot_areaerrorbar(real(solind_C(cohort1_id,:)),options)
hold on
options.color_area = [59 164 209]/255;
options.color_line = [59 164 209]/255;
plot_areaerrorbar(real(solind_C(cohort2_id,:)),options)
xlabel('Time')
ylabel('Cancer cells C, cells')
set(gca,'FontSize',18)
ylim([0 8e6])

% plot mean and std for PAC-1 for each sub cohort
fig1 = figure
hold on
options.handle = fig1;
options.color_area = [209,105,61]/255;
options.color_line = [209,105,61]/255;
options.alpha = 0.3;
options.error = 'std';
options.line_width = 2;
options.x_axis = time;
plot_areaerrorbar(solind_PAC1(cohort1_id,:),options)
hold on
options.color_area = [59 164 209]/255;
options.color_line = [59 164 209]/255;
plot_areaerrorbar(solind_PAC1(cohort2_id,:),options)
xlabel('Time')
ylabel('PAC-1')
set(gca,'FontSize',18)
ylim([0 2500])

% plot mean and std for TRAIL for each sub cohort
fig1 = figure
options.handle = fig1;
options.color_area = [209,105,61]/255;
options.color_line = [209,105,61]/255;
options.alpha = 0.3;
options.error = 'std';
options.line_width = 2;
options.x_axis = time;
plot_areaerrorbar(solind_TRAIL(cohort1_id,:),options)
hold on
options.color_area = [59 164 209]/255;
options.color_line = [59 164 209]/255;
plot_areaerrorbar(solind_TRAIL(cohort2_id,:),options)
xlabel('Time')
ylabel('TRAIL')
set(gca,'FontSize',18)
ylim([0 400])

%% load optimal protocols and plot cohort dynamics for High PAC-1 elimination cohort and 50 virtual patients from Normal PAC-1 elimination cohort
load('Results_Nov8')

OptimalDose = Results{end}.OptimalDose;   

%setting colour map for the heatmap
sumX = [247,252,240
224,243,219
204,235,197
168,221,181
123,204,196
78,179,211
43,140,190
8,104,172
8,64,129]/255;

%plotting optimal PAC-1 treatment for patients in the patients in the High PAC-1
%elimination sub cohort
figure
h = heatmap(OptimalDose(cohort1_id,1:11),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')
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
h.YDisplayLabels(1:length(cohort1_id)) = {""}

s = struct(h);
s.XAxis.TickLabelRotation = 0;  % vertical
title('PAC-1')

%plotting optimal PAC-1 treatment for the first 50 patients in the patients in the Normal PAC-1
%elimination sub cohort
figure
h = heatmap(OptimalDose(cohort2_id(1:50),1:11),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')
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
h.YDisplayLabels(1:50) = {""}

s = struct(h);
s.XAxis.TickLabelRotation = 0;  % vertical
title('Cohort 2')

%plotting optimal TRAIL treatment for patients in the patients in the High PAC-1
%elimination sub cohort
figure
h = heatmap(OptimalDose(cohort1_id,12:22),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')
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

h.YDisplayLabels(1:length(cohort1_id)) = {""}

s = struct(h);
s.XAxis.TickLabelRotation = 0;  % vertical
title('PAC-1')

%plotting optimal TRAIL treatment for the first 50 patients in the patients in the Normal PAC-1
%elimination sub cohort
figure
h = heatmap(OptimalDose(cohort2_id(1:50),12:22),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')
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
h.YDisplayLabels(1:50) = {""}

s = struct(h);
s.XAxis.TickLabelRotation = 0;  % vertical
title('Cohort 2')
