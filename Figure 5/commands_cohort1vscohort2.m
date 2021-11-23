% simulate average model behaviour

load VP.mat
load VP_May31.mat
VirtualP=[VP; VP_May31];

load('cohort1_id.mat')
cohort2_id = [];
for ii = 1:400
    if isempty(find(cohort1_id == ii))==1
        cohort2_id = [cohort2_id;ii];  
    end
end

cohort1_parameters = VirtualP(cohort1_id,:);
cohort2_parameters = VirtualP(cohort2_id,:);


% plotting cohort parameters
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

% plotting cohort parameters
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

% plotting cohort parameters
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

% plotting cohort parameters
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

% plotting cohort parameters
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

% plotting cohort parameters
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

% plotting cohort parameters
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
    solind_C(i,:) = deval(sol,time,1);
    solind_D(i,:) = deval(sol,time,2);
    solind_PAC1(i,:) = deval(sol,time,3);
    solind_TRAIL(i,:) = deval(sol,time,4);
    
    CumulativeTumourBurden_mat(i,:) = trapz(time,solind_C(i,:));
    
end


figure
histogram(log(real(CumulativeTumourBurden_mat(cohort1_id,:))),10,'FaceColor',[209,105,61]/255,'EdgeColor','none')
hold on
histogram(log(real(CumulativeTumourBurden_mat(cohort2_id,:))),100,'FaceColor',[59 164 209]/255,'EdgeColor','none')
xlabel('log(Cumulative Tumour Burden)')
ylabel('Frequency') 
set(gca,'FontSize',18)
box off
legend('Cohort 1','Cohort 2')


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

fig2 = figure
hold on
options.handle = fig2;
options.color_area = [209,105,61]/255;
options.color_line = [209,105,61]/255;
options.alpha = 0.3;
options.error = 'std';
options.line_width = 2;
options.x_axis = time;
plot_areaerrorbar(real(solind_D(cohort1_id,:)),options)
hold on
options.color_area = [59 164 209]/255;
options.color_line = [59 164 209]/255;
plot_areaerrorbar(real(solind_D(cohort2_id,:)),options)
xlabel('Time')
ylabel('Dead cells D, cells')
set(gca,'FontSize',18)

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

%% load optimal protocols and plot cohort dynamics
load('Results_Nov8')

OptimalDose = Results{end}.OptimalDose;   


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

figure
h = heatmap(OptimalDose(cohort2_id,1:11),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')
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
h.YDisplayLabels(1:length(cohort2_id)) = {""}

s = struct(h);
s.XAxis.TickLabelRotation = 0;  % vertical
title('Cohort 2')

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

figure
h = heatmap(OptimalDose(cohort2_id,12:22),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')
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
h.YDisplayLabels(1:length(cohort2_id)) = {""}

s = struct(h);
s.XAxis.TickLabelRotation = 0;  % vertical
title('Cohort 2')
