% determine a sub cohort level optimal protocol, i.e. a single protocol
% that is optimal for all individuals in the High PAC-1 elimination sub
% group and then a protocol that is optimal for all individuals in the
% Normal PAC-1 elimination sub group

% load virtual patient parameters
load newcohorts.mat

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
OptimalDose = zeros(size(200,1),nvars); %Matrix to store the optimal doses for each patient


%% optimising protocol for each sub cohort

% genetic algorithm for optimising protocol for High PAC-1 elimination
% patients
FObjective=@(x)ProtocolOptimizerObjectiveCOHORTS_V2(x,p,newpatients_C1);
opts = optimoptions('ga','EliteCount',15,'MaxStallGenerations',4,'FunctionTolerance',1e-6, 'MaxGenerations',100,'Display','iter','PopulationSize',60);
[x,fval_C1] = ga(FObjective, nvars,[],[],[],[],LB,UB,[],IntCon,opts);   

%Store the optimal dosing regime ****
OptimalDose_C1 = x;
Results.OptimalDose_C1=OptimalDose_C1;
Results.fval=fval_C1;

% genetic algorithm for optimising protocol for Normal PAC-1 elimination
% patients
FObjective=@(x)ProtocolOptimizerObjectiveCOHORTS_V2(x,p,newpatients_C2);
opts = optimoptions('ga','EliteCount',15,'MaxStallGenerations',4,'FunctionTolerance',1e-6, 'MaxGenerations',100,'Display','iter','PopulationSize',60);
[x,fval_C2] = ga(FObjective, nvars,[],[],[],[],LB,UB,[],IntCon,opts);   
  
%Store the optimal dosing regime ****
OptimalDose_C2 = x;
Results.OptimalDose_C2=OptimalDose_C2;
Results.fval=fval_C2;

%saving optimal doses
save('Results_NewSubcohortOptimal2','OptimalDose_C2','OptimalDose_C1')

%% plotting optimal sub cohort protocols

PAC1COL = [255,245,235
254,230,206
253,208,162
253,174,107
253,141,60
241,105,19
217,72,1
166,54,3
127,39,4]/255;

TRAILCOL = [252,251,253
239,237,245
218,218,235
188,189,220
158,154,200
128,125,186
106,81,163
84,39,143
63,0,125]/255;

colgrid = linspace(0,10.1,8)

figure
hold on 
for i = 1:11
    B = bar(i-1,OptimalDose_C1(i))
    colindex = find(OptimalDose_C1(i)<=colgrid,1);
    set(B,'FaceColor',PAC1COL(colindex,:),'EdgeColor','none')
end
c=colorbar
colormap(PAC1COL)
xlabel('Dosage week')
ylabel('Dosage increase')
title('PAC-1')
set(gca,'FontSize',18)
ylim([0 10])
set(c,'ticks',[0 0.5 1],'tickLabels',{'0','5','10'})
    
figure
hold on 
for i = 12:22
    B = bar(i-12,OptimalDose_C1(i))
    colindex = find(OptimalDose_C1(i)<=colgrid,1);
    set(B,'FaceColor',TRAILCOL(colindex,:),'EdgeColor','none')
end
c=colorbar
colormap(TRAILCOL)
xlabel('Dosage week')
ylabel('Dosage increase')
set(gca,'FontSize',18)
title('TRAIL')
ylim([0 10])
set(c,'ticks',[0 0.5 1],'tickLabels',{'0','5','10'})

figure
hold on 
for i = 1:11
    B = bar(i-1,OptimalDose_C2(i))
    colindex = find(OptimalDose_C2(i)<=colgrid,1);
    set(B,'FaceColor',PAC1COL(colindex,:),'EdgeColor','none')
end
c=colorbar
colormap(PAC1COL)
xlabel('Dosage week')
ylabel('Dosage increase')
title('PAC-1')
set(gca,'FontSize',18)
ylim([0 10])
set(c,'ticks',[0 0.5 1],'tickLabels',{'0','5','10'})
    
figure
hold on 
for i = 12:22
    B = bar(i-12,OptimalDose_C2(i))
    colindex = find(OptimalDose_C2(i)<=colgrid,1);
    set(B,'FaceColor',TRAILCOL(colindex,:),'EdgeColor','none')
end
c = colorbar
colormap(TRAILCOL)
xlabel('Dosage week')
ylabel('Dosage increase')
set(gca,'FontSize',18)
title('TRAIL')
ylim([0 10])
set(c,'ticks',[0 0.5 1],'tickLabels',{'0','5','10'})

%% simulate individuals in each cohort under optimal protocol and plot dynamics


for j=1:200 %For each patient, find the optimal dosing regime
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
    
    %simulate model on optimal dosage and store
    [sol] = modelsimulator_multidose(OptimalDose_C1,p);   
    SolMat_C1{j} = sol;
        
end
    save('Results_C1','Results','SolMat_C1')


for j=201:400 %For each patient, find the optimal dosing regime
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
    
    %simulate model on optimal dosage and store
    [sol] = modelsimulator_multidose(OptimalDose_C2,p);   
    SolMat_C2{j} = sol;
        
end

    save('Results_C2','Results','SolMat_C2')

%% plotting comparison between subcohort 1 and subcohort 2 under optimal protocols

% plotting comparison between mean and std of cohort 1 and cohort 2 for
% number of cells


for i = 1:200
    
    sol = SolMat_C1{i};
    time = linspace(0,p.tf,1000);
    solind_C_C1(i,:) = deval(sol,time,1);
    solind_D_C1(i,:) = deval(sol,time,2);
    solind_PAC1_C1(i,:) = deval(sol,time,3);
    solind_TRAIL_C1(i,:) = deval(sol,time,4);
    
    CumulativeTumourBurden_mat_cohort1(i,:) = trapz(time,solind_C_C1(i,:));
    
end


for i = 201:400
    
    sol = SolMat_C2{i};
    time = linspace(0,p.tf,1000);
    solind_C_C2(i,:) = deval(sol,time,1);
    solind_D_C2(i,:) = deval(sol,time,2);
    solind_PAC1_C2(i,:) = deval(sol,time,3);
    solind_TRAIL_C2(i,:) = deval(sol,time,4);
    
    CumulativeTumourBurden_mat_cohort2(i,:) = trapz(time,solind_C_C2(i,:));
    
end


figure
histogram(log(real(CumulativeTumourBurden_mat_cohort1)),200,'FaceColor',[47 85 151]/255,'EdgeColor','none')
hold on
histogram(log(real(CumulativeTumourBurden_mat_cohort2)),200,'FaceColor',[84 130 53]/255,'EdgeColor','none')
xlabel('log(Cumulative Tumour Burden)')
ylabel('Frequency') 
set(gca,'FontSize',18)
box off


fig1 = figure
hold on
plot(time,solind_C_C1,'Color',[47 85 151]/255)
plot(time,solind_C_C2,'Color',[84 130 53]/255)
xlabel('Time')
ylabel('Cancer cells C, log(cells)')
set(gca,'FontSize',18)
set(gca,'yscale','log')

fig1 = figure
hold on
options.handle = fig1;
options.color_area = [47 85 151]/255;
options.color_line =  [47 85 151]/255;
options.alpha = 0.3;
options.error = 'std';
options.line_width = 2;
options.x_axis = time;
plot_areaerrorbar(real(solind_C_C1),options)
hold on
options.color_area = [84 130 53]/255;
options.color_line =  [84 130 53]/255;
plot_areaerrorbar(real(solind_C_C2),options)
xlabel('Time')
ylabel('Cancer cells C, cells')
set(gca,'FontSize',18)
ylim([0 2e6])
hold on 
l1 = plot(time, mean(real(solind_C_C1)),'LineWidth',2,'Color',[47 85 151]/255)
hold on 
l2 = plot(time, mean(real(solind_C_C2)),'LineWidth',2,'Color',[84 130 53]/255)
legend([l1 l2],{'Cohort 1','Cohort 2'})


fig1 = figure
hold on
options.handle = fig1;
options.color_area = [47 85 151]/255;
options.color_line =  [47 85 151]/255;
options.alpha = 0.3;
options.error = 'std';
options.line_width = 2;
options.x_axis = time;
plot_areaerrorbar(log(real(solind_C_C1)),options)
hold on 
options.color_area = [84 130 53]/255;
options.color_line =  [84 130 53]/255;
plot_areaerrorbar(log(real(solind_C_C2)),options)
xlabel('Time')
ylabel('Cancer cells C, log_{10}(cells)')
set(gca,'FontSize',18)


figure
hold on 
errorbar(time,mean(log(real(solind_C_C1))),std(log(real(solind_C_C1))))
errorbar(time,mean(log(real(solind_C_C2))),std(log(real(solind_C_C2))))

%% plotting a heatmap of patient tumour volume over time coloured by cells <0 and otherwise
load('Results_C2')

sumX = [202 240 248
    173 232 244
     144 224 239
    72 202 220
    6 180 216
    0 150 199
    0 119 182
    2 62 138
    3 4 94]/255;


for i = 1:200
    
    sol = SolMat_C1{i};
    times = linspace(0,p.tf,100);
    solind_C_C1(i,:) = deval(sol,times,1);
end

Csol_C1 = real(solind_C_C1);
for pcount = 1:200
    locs = find(Csol_C1(pcount,:)<1e-1);
    Csol_C1(pcount,locs) = 0;
end

colvec = [10 1e2 1e3 1e4 1e5 1e6 1e7+1];

[A,B] = sort(Csol_C1(:,end))

figure
hold on
for ii = 1:200
    for tcount = 1:100
      if Csol_C1(ii,tcount) >0 && Csol_C1(ii,tcount)<1e7
          locval=find(colvec>=Csol_C1(ii,tcount),1);
          colval = sumX(locval,:);
      elseif Csol_C1(ii,tcount) <= 0
          colval = sumX(1,:);
      elseif Csol_C1(ii,tcount)>=1e7
          colval = sumX(end,:);
      end
      fill([tcount-1 tcount-1 tcount tcount],[ii-1 ii ii ii-1],colval,'EdgeColor','none');
        
%h = heatmap(Csol_C1(:,1:100),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')

    end
end

c= colorbar
colormap(sumX)
ylim([0 200])
xlabel('Time (days)')
ylabel('High PAC-1 elimination sub cohort')
set(gca,'FontSize',25)
set(gca,'xtick',[0 50 100],'xticklabels',{'0','35','70'})
set(c,'ticks',[0 0.1111 0.2222 0.3333 0.4444 0.5555 0.6666 0.777 0.8889 1],'tickLabels',{'0','1','10','10^2','10^3','10^4','10^5','10^6','10^7','10^8'})
set(gca,'yticklabels',[''])

sumX = [216 243 220
    183 228 199
    149 213 178
    116 198 157
    81 183 135
    64 145 108
    45 106 79
    27 67 50
    8 28 21]/255;


for i = 201:400
    
    sol = SolMat_C2{i};
    times = linspace(0,p.tf,100);
    solind_C_C2(i-200,:) = deval(sol,times,1);
end

Csol_C2= real(solind_C_C2);
for pcount = 1:200
    locs = find(Csol_C2(pcount,:)<1e-1);
    Csol_C2(pcount,locs) = 0;
end

colvec = [10 1e2 1e3 1e4 1e5 1e6 1e7+1];

figure
hold on
for ii = 1:200
    for tcount = 1:100
      if Csol_C2(ii,tcount) >0 && Csol_C2(ii,tcount)<1e7
          locval=find(colvec>=Csol_C2(ii,tcount),1);
          colval = sumX(locval,:);
      elseif Csol_C2(ii,tcount) <= 0
          colval = sumX(1,:);
      elseif Csol_C2(ii,tcount)>=1e7
          colval = sumX(end,:);
      end
      fill([tcount-1 tcount-1 tcount tcount],[ii-1 ii ii ii-1],colval,'EdgeColor','none');
        
%h = heatmap(Csol_C1(:,1:100),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')

    end
end

c= colorbar
colormap(sumX)
ylim([0 200])
xlabel('Time (days)')
ylabel('Normal PAC-1 elimination sub cohort')
set(gca,'FontSize',25)
set(gca,'xtick',[0 50 100],'xticklabels',{'0','35','70'})
set(c,'ticks',[0 0.1111 0.2222 0.3333 0.4444 0.5555 0.6666 0.777 0.8889 1],'tickLabels',{'0','1','10','10^2','10^3','10^4','10^5','10^6','10^7','10^8'})
set(gca,'yticklabels',[''])

%% plotting a heatmap of patient tumour volume over time coloured by cells <0 and otherwise
load('Results_C2')

sumX = [202 240 248
    173 232 244
     144 224 239
    72 202 220
    6 180 216
    0 150 199
    0 119 182
    2 62 138
    3 4 94]/255;


for i = 1:200
    
    sol = SolMat_C1{i};
    times = linspace(0,p.tf,100);
    solind_C_C1(i,:) = deval(sol,times,1);
end

Csol_C1 = real(solind_C_C1);
for pcount = 1:200
    locs = find(Csol_C1(pcount,:)<1e-1);
    Csol_C1(pcount,locs) = 0;
end

colvec = [10 1e2 1e3 1e4 1e5 1e6 1e7+1];


%calculate AUC for each patient
AUC = sum(real(solind_C_C1)')*(times(2)-times(1));
[A,B] = sort(AUC);

figure
hold on
for ii = 1:200
    for tcount = 1:100
      if Csol_C1(B(ii),tcount) >0 && Csol_C1(B(ii),tcount)<1e7
          locval=find(colvec>=Csol_C1(B(ii),tcount),1);
          colval = sumX(locval,:);
      elseif Csol_C1(B(ii),tcount) <= 0
          colval = sumX(1,:);
      elseif Csol_C1(B(ii),tcount)>=1e7
          colval = sumX(end,:);
      end
      fill([tcount-1 tcount-1 tcount tcount],[ii-1 ii ii ii-1],colval,'EdgeColor','none');
        
%h = heatmap(Csol_C1(:,1:100),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')

    end
end

c= colorbar
colormap(sumX)
ylim([0 200])
xlabel('Time (days)')
ylabel('High PAC-1 elimination sub cohort')
set(gca,'FontSize',18)
set(gca,'xtick',[0 50 100],'xticklabels',{'0','35','70'})
set(c,'ticks',[0 0.1111 0.2222 0.3333 0.4444 0.5555 0.6666 0.777 0.8889 1],'tickLabels',{'0','1','10','10^2','10^3','10^4','10^5','10^6','10^7','10^8'})


sumX = [216 243 220
    183 228 199
    149 213 178
    116 198 157
    81 183 135
    64 145 108
    45 106 79
    27 67 50
    8 28 21]/255;


for i = 201:400
    
    sol = SolMat_C2{i};
    times = linspace(0,p.tf,100);
    solind_C_C2(i-200,:) = deval(sol,times,1);
end

Csol_C2= real(solind_C_C2);
for pcount = 1:200
    locs = find(Csol_C2(pcount,:)<1e-1);
    Csol_C2(pcount,locs) = 0;
end

colvec = [10 1e2 1e3 1e4 1e5 1e6 1e7+1];

%calculate AUC for each patient
AUC = sum(real(solind_C_C2)')*(times(2)-times(1));
[A,B] = sort(AUC);

figure
hold on
for ii = 1:200
    for tcount = 1:100
      if Csol_C2(B(ii),tcount) >0 && Csol_C2(B(ii),tcount)<1e7
          locval=find(colvec>=Csol_C2(B(ii),tcount),1);
          colval = sumX(locval,:);
      elseif Csol_C2(B(ii),tcount) <= 0
          colval = sumX(1,:);
      elseif Csol_C2(B(ii),tcount)>=1e7
          colval = sumX(end,:);
      end
      fill([tcount-1 tcount-1 tcount tcount],[ii-1 ii ii ii-1],colval,'EdgeColor','none');
        
%h = heatmap(Csol_C1(:,1:100),'CellLabelColor','none','Colormap',sumX,'GridVisible','off')

    end
end

c= colorbar
colormap(sumX)
ylim([0 200])
xlabel('Time (days)')
ylabel('Normal PAC-1 elimination sub cohort')
set(gca,'FontSize',18)
set(gca,'xtick',[0 50 100],'xticklabels',{'0','35','70'})
set(c,'ticks',[0 0.1111 0.2222 0.3333 0.4444 0.5555 0.6666 0.777 0.8889 1],'tickLabels',{'0','1','10','10^2','10^3','10^4','10^5','10^6','10^7','10^8'})

%%
