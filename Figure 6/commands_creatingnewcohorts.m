%% code to generate a new cohort of virtual patients based on the sub cohorts identified

%load virtual patient parameters
load VP.mat
load VP_May31.mat
VirtualP=[VP; VP_May31];

%load subcohort ides
load('cohort1_id.mat')
cohort2_id = [];
for ii = 1:400
    if isempty(find(cohort1_id == ii))==1
        cohort2_id = [cohort2_id;ii];  
    end
end

cohort1_parameters = VirtualP(cohort1_id,:);
cohort2_parameters = VirtualP(cohort2_id,:);

%% estimating PDF for Cohort 1 (High PAC-1 elimination cohort)

%estimating and resampling from the PDF for cohort 1 for parameter Vd_PAC1
pd_C1_P1 = fitdist(cohort1_parameters(:,1),'kernel');
Y = abs(random(pd_C1_P1, [1000,1]));
figure
histogram(cohort1_parameters(:,1),'Normalization','probability')
hold on 
histogram(Y,'Normalization','probability')

%estimating and resampling from the PDF for cohort 1 for parameter CL_PAC1
pd_C1_P2 = fitdist(cohort1_parameters(:,2),'kernel');
Y = abs(random(pd_C1_P2, [1000,1]));
figure
histogram(cohort1_parameters(:,2),'Normalization','probability')
hold on 
histogram(Y,'Normalization','probability')

%estimating and resampling from the PDF for cohort 1 for parameter Vd_TRAIL
pd_C1_P3 = fitdist(cohort1_parameters(:,3),'kernel');
Y = abs(random(pd_C1_P3, [1000,1]));
figure
histogram(cohort1_parameters(:,3),'Normalization','probability')
hold on 
histogram(Y,'Normalization','probability')

%estimating and resampling from the PDF for cohort 1 for parameter CL_TRAIL
pd_C1_P4 = fitdist(cohort1_parameters(:,4),'kernel');
Y = abs(random(pd_C1_P4, [1000,1]));
figure
histogram(cohort1_parameters(:,4),'Normalization','probability')
hold on 
histogram(Y,'Normalization','probability')

%estimating and resampling from the PDF for cohort 1 for parameter BW
pd_C1_P5 = fitdist(cohort1_parameters(:,5),'kernel');
Y = abs(random(pd_C1_P5, [1000,1]));
figure
histogram(cohort1_parameters(:,5),'Normalization','probability')
hold on 
histogram(Y,'Normalization','probability')

%% estimating PDF for cohort 2 (Normal PAC-1 elimination cohort)

%estimating and resampling from the PDF for cohort 2 for parameter Vd_PAC1
pd_C2_P1 = fitdist(cohort2_parameters(:,1),'kernel');
Y = random(pd_C2_P1, [1000,1]);
figure
histogram(cohort2_parameters(:,1),'Normalization','probability')
hold on 
histogram(Y,'Normalization','probability')

pd_C2_P2 = fitdist(cohort2_parameters(:,2),'kernel');
Y = random(pd_C2_P2, [1000,1]);
figure
histogram(cohort2_parameters(:,2),'Normalization','probability')
hold on 
histogram(Y,'Normalization','probability')

pd_C2_P3 = fitdist(cohort2_parameters(:,3),'kernel');
Y = random(pd_C2_P3, [1000,1]);
figure
histogram(cohort2_parameters(:,3),'Normalization','probability')
hold on 
histogram(Y,'Normalization','probability')

pd_C2_P4 = fitdist(cohort2_parameters(:,4),'kernel');
Y = random(pd_C2_P4, [1000,1]);
figure
histogram(cohort2_parameters(:,4),'Normalization','probability')
hold on 
histogram(Y,'Normalization','probability')

pd_C2_P5 = fitdist(cohort2_parameters(:,5),'kernel');
Y = random(pd_C2_P5, [1000,1]);
figure
histogram(cohort2_parameters(:,5),'Normalization','probability')
hold on 
histogram(Y,'Normalization','probability')

%% resample patients in each new cohort
   
    newpatients_C1 = abs([random(pd_C1_P1,  [200,1]),random(pd_C1_P2,  [200,1]),random(pd_C1_P3,  [200,1]),random(pd_C1_P4, [200,1]),random(pd_C1_P5,  [200,1])]);
    newpatients_C2 = abs([random(pd_C2_P1,  [200,1]),random(pd_C2_P2,  [200,1]),random(pd_C2_P3,  [200,1]),random(pd_C2_P4,  [200,1]),random(pd_C2_P5, [200,1])]);
    
%% Plot box plots of the parameters to show how they relate

figure
subplot(2,5,1)
notBoxPlot([[cohort1_parameters(:,1);NaN(200-35,1)],newpatients_C1(:,1)])
set(gca,'FontSize',18)
set(gca,'xtick',[1,2],'xticklabels',{'Original','New'})
ylabel('VD_{PAC1}')
title('PAC-1 volume distribution')

subplot(2,5,2)
notBoxPlot([[cohort1_parameters(:,2);NaN(200-35,1)],newpatients_C1(:,2)])
set(gca,'FontSize',18)
set(gca,'xtick',[1,2],'xticklabels',{'Original','New'})
ylabel('CL_{PAC1}')
title('PAC-1 clearance')

subplot(2,5,3)
notBoxPlot([[cohort1_parameters(:,3);NaN(200-35,1)],newpatients_C1(:,3)])
set(gca,'FontSize',18)
set(gca,'xtick',[1,2],'xticklabels',{'Original','New'})
ylabel('VD_{TRAIL}')
title('TRAIL Volume distribution')

subplot(2,5,4)
notBoxPlot([[cohort1_parameters(:,4);NaN(200-35,1)],newpatients_C1(:,4)])
set(gca,'FontSize',18)
set(gca,'xtick',[1,2],'xticklabels',{'Original','New'})
ylabel('CL_{TRAIL}')
title('TRAIL clearance')

subplot(2,5,5)
notBoxPlot([[cohort1_parameters(:,5);NaN(200-35,1)],newpatients_C1(:,5)])
set(gca,'FontSize',18)
set(gca,'xtick',[1,2],'xticklabels',{'Original','New'})
ylabel('BW (kg)')
title('Body weight')

subplot(2,5,6)
notBoxPlot([cohort2_parameters(:,1),[newpatients_C2(:,1);NaN(165,1)]])
set(gca,'FontSize',18)
set(gca,'xtick',[1,2],'xticklabels',{'Original','New'})
ylabel('VD_{PAC1}')
title('PAC-1 volume distribution')

subplot(2,5,7)
notBoxPlot([cohort2_parameters(:,2),[newpatients_C2(:,2);NaN(165,1)]])
set(gca,'FontSize',18)
set(gca,'xtick',[1,2],'xticklabels',{'Original','New'})
ylabel('CL_{PAC1}')
title('PAC-1 clearance')

subplot(2,5,8)
notBoxPlot([cohort2_parameters(:,3),[newpatients_C2(:,3);NaN(165,1)]])
set(gca,'FontSize',18)
set(gca,'xtick',[1,2],'xticklabels',{'Original','New'})
ylabel('VD_{TRAIL}')
title('TRAIL Volume distribution')

subplot(2,5,9)
notBoxPlot([cohort2_parameters(:,4),[newpatients_C2(:,4);NaN(165,1)]])
set(gca,'FontSize',18)
set(gca,'xtick',[1,2],'xticklabels',{'Original','New'})
ylabel('CL_{TRAIL}')
title('TRAIL clearance')

subplot(2,5,10)
notBoxPlot([cohort2_parameters(:,5),[newpatients_C2(:,5);NaN(165,1)]])
set(gca,'FontSize',18)
set(gca,'xtick',[1,2],'xticklabels',{'Original','New'})
ylabel('BW (kg)')
title('Body weight')


%%

% plotting cohort parameters
figure 
v = violinplot([newpatients_C1(:,1),newpatients_C2(:,1)]);   
ylabel('VD_{PAC1}')
title('PAC-1 volume distribution')
set(gca,'xtick',[1,2],'xticklabels',{'New cohort 1','New cohort 2'})
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;

figure 
v = violinplot([newpatients_C1(:,2),newpatients_C2(:,2)]);   
ylabel('CL_{PAC1}')
title('PAC-1 clearance')
set(gca,'xtick',[1,2],'xticklabels',{'New cohort 1','New cohort 2'})
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;


figure 
v = violinplot([newpatients_C1(:,3),newpatients_C2(:,3)]);   
ylabel('VD_{TRAIL}') 
title('TRAIL volume distribution')
set(gca,'xtick',[1,2],'xticklabels',{'New cohort 1','New cohort 2'})
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;


figure 
v = violinplot([newpatients_C1(:,4),newpatients_C2(:,4)]);   
ylabel('CL_{TRAIL}')
title('TRAIL clearance')
set(gca,'xtick',[1,2],'xticklabels',{'New cohort 1','New cohort 2'})
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;

figure 
v = violinplot([newpatients_C1(:,5),newpatients_C2(:,5)]);   
ylabel('BW (kg)')
title('Body weight')
set(gca,'xtick',[1,2],'xticklabels',{'New cohort 1','New cohort 2'})
set(gca,'FontSize',18);
v(1).ViolinColor = [209,105,61]/255;
v(2).ViolinColor = [59 164 209]/255;
v(1).ViolinAlpha = 0.5;
v(2).ViolinAlpha = 0.5;
%%

%save('newcohorts.mat','newpatients_C1','newpatients_C2')
