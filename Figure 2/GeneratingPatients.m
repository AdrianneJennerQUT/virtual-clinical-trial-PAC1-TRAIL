function VP=GeneratingPatients(num_patients)

% Calculate BW
load PopulationDataWomen
load PopulationDataMen
possibleGenders={'Man','Woman'};
possibleCountries=PopulationDataWomen.Properties.RowNames; %get a list of all countries
%for which we have BMI and height data, this could be restricted to a
%single country or a single region as well, use the ISO in the
%PopulationDataWomen or PopulationDataMen files as country names.

np=num_patients; %number of patients to simulate
% nGen=randi([1, 2],1,np); %construct random sequence of 1 and 2 to construct 
%sequence of gender information about population
% GenderPop=possibleGenders(nGen); %generate distribution of genders within population
GenderPop=cell(1,np);
%GenderPop(:)=cellstr('Man');
GenderPop(:)=cellstr('Woman') ;
%Construct population from all possible countries
randCount=randi(length(possibleCountries),1,np); %construct random sequence from 1 to total 
%number of countries listed 
Countries=possibleCountries(randCount);%generate distribution of countries for each individual

% %Construct population from one single country
 Countries=cell(1,np);
 Countries(:)=cellstr('USA');

MeanBMI=zeros(1,np);
varBMI=zeros(1,np);
MeanHeight=zeros(1,np);
varHeight=zeros(1,np);

for nn=1:np
    if strcmp(GenderPop(nn),'Man')==1
        MeanBMI(nn)=table2array(PopulationDataMen(Countries(nn),'MeanBMI'));
        varBMI(nn)=table2array(PopulationDataMen(Countries(nn),'CalculatedSDforBMI'));
        MeanHeight(nn)=table2array(PopulationDataMen(Countries(nn),'Meanheightcm'));
        varHeight(nn)=table2array(PopulationDataMen(Countries(nn),'CalculatedSDforHeight'));
    else
        MeanBMI(nn)=table2array(PopulationDataWomen(Countries(nn),'MeanBMI'));
        varBMI(nn)=table2array(PopulationDataWomen(Countries(nn),'CalculatedSDforBMI'));
        MeanHeight(nn)=table2array(PopulationDataWomen(Countries(nn),'Meanheightcm'));
        varHeight(nn)=table2array(PopulationDataWomen(Countries(nn),'CalculatedSDforHeight'));
    end
end

%Random sampling of BMIs and heights according to country and gender
%constructed from a normal distribution
BMIs=normrnd(MeanBMI,varBMI);%kg/m^2
Heights=normrnd(MeanHeight,varHeight)/100;%m
     
%Physiological characteristics of population
BW=BMIs.*Heights.^2; %average body weight in kgs

figure
histogram(BMIs,'FaceColor',[6 141 157]/255,'EdgeColor',[1 1 1])
ylabel('Frequency')
xlabel('BMI')
grid on
set(gca,'FontSize',18)


figure
histogram(Heights,'FaceColor',[128 222 217]/255,'EdgeColor',[1 1 1])
ylabel('Frequency')
xlabel('H')
grid on
set(gca,'FontSize',18)


figure
histogram(BW,'FaceColor',[174 236 239]/255,'EdgeColor',[1 1 1])
ylabel('Frequency')
xlabel('BW')
grid on 
set(gca,'FontSize',18)

% PAC-1 parameters
%Fixed effects
theta_VdcPAC1=0.7;%Vd central (L/kg)
%theta_VdePAC1=2.7;%Vd terminal (L/kg)
%theta_VdssPAC1=1.6;%Vd steady state (L/kg)
%theta_thalf_alphaPAC1=0.4;%distribution half-life (h)
%theta_thalf_betaPAC1=3.1;%terminal half-life (h)
theta_ClPAC1=0.5744;%total body clearance (L/kg/h)

%Standard deviations
SD_VdcPAC1=0.3;
%SD_VdePAC1=0.9;
%SD_VdssPAC1=0.5;
%SD_thalf_alphaPAC1=0.2;
%SD_thalf_betaPAC1=0.7;
SD_ClPAC1=0.0757;

% FixedValuesPAC1=[theta_VdcPAC1 theta_VdePAC1 theta_VdssPAC1...
%     theta_thalf_alphaPAC1 theta_thalf_betaPAC1 theta_ClPAC1];
% VarPAC1 = [SD_VdcPAC1^2 0 0 0 0 0;
%            0 SD_VdePAC1^2 0 0 0 0;
%            0 0 SD_VdssPAC1^2 0 0 0;
%            0 0 0 SD_thalf_alphaPAC1^2 0 0;
%            0 0 0 0 SD_thalf_betaPAC1^2 0;
%            0 0 0 0 0 SD_ClPAC1^2];

FixedValuesPAC1=[theta_VdcPAC1 theta_ClPAC1];
VarPAC1 = [SD_VdcPAC1^2 0
           0 SD_ClPAC1^2];

%Sample patients
rng(1);%to make sure you generate the same patients each time you run this (remove to make new cohort) - seed
PatientParamsPAC1=mvnrnd(FixedValuesPAC1,VarPAC1,num_patients); %randomly generating patients

% TRAIL parameters
%Fixed effects
theta_Vd_TRAIL=4.28;%Vd (L)
theta_Cl_TRAIL=116;%clearance (L/day)

%Standard deviations (assuming % is coefficient of variation)
SD_Vd_TRAIL=0.279*theta_Vd_TRAIL;
SD_Cl_TRAIL=0.326*theta_Cl_TRAIL;

FixedValuesTRAIL=[theta_Vd_TRAIL theta_Cl_TRAIL];
VarTRAIL=[SD_Vd_TRAIL^2 0;
          0 SD_Cl_TRAIL];

%Sample patients
PatientParamsTRAIL=mvnrnd(FixedValuesTRAIL,VarTRAIL,num_patients); %randomly generating patients

figure
s = scatterhist(PatientParamsPAC1(:,1),PatientParamsPAC1(:,2),'Kernel','off','Location','SouthEast','Direction','out','Color',[83 89 154]/255,'LineWidth',[15])
set(gca,'FontSize',18)
xlabel('VD_{PAC1}')
ylabel('CL_{PAC1}')
figure
s = scatterhist(PatientParamsTRAIL(:,1),PatientParamsTRAIL(:,2),'Kernel','off','Location','SouthEast','Direction','out','Color',[109 157 197]/255,'LineWidth',[15])
set(gca,'FontSize',18)
xlabel('VD_{TRAIL}')
ylabel('CL_{TRAIL}')


VP=[PatientParamsPAC1,PatientParamsTRAIL,BW'];
