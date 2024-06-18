addpath('../Functions/');
load('../Data/PWHPrevData.mat');
load('../Data/PWHNewDiagData.mat');
load('../Data/lifeTable.mat');

%read data
firstYearData=PWHPrevData(PWHPrevData.Year==2008,:);
ages=[0;firstYearData.AgeGroup;100];

%set some parameter values
minAge=1;
maxAge=101;

%death rate
mu=-log(1-.025);

%this is the max likelihood fitting to the 2008 data to initialize the
%population
[maxMu,maxSigma,popDist,~]=...
    age_Distribution_MaxLikelihood(ages,[0;firstYearData.Cases]);

%time step
dt=0.25;

%age mesh initialization
ageMesh=linspace(minAge,maxAge,(maxAge-minAge)/dt)';

%system matrices
A=assemble_Age_Matrix(minAge,maxAge,dt);
M=assemble_Mass_Matrix(size(A,1),dt);
Mu=assemble_Mortality_Matrix(lifeTable,ageMesh);

%initial distribution
P_0=sum(firstYearData.Cases)*popDist(ageMesh);

sols=[P_0];

nTimeSteps=16/dt;
PLast=P_0;

%beginningYear
yearStart=2008;
timeOutput=[yearStart];

curTime=yearStart+dt;
timeOutput=[timeOutput;curTime];
curYearData=PWHNewDiagData(PWHNewDiagData.Year==yearStart,:);  
[diagMu,diagSigma,newDiags,~]=...
    age_Distribution_MaxLikelihood(ages,[0;curYearData.Cases]);
newEntries=sum(curYearData.Cases)*newDiags(ageMesh);

%initialize with Heun's method
PInt=PLast+dt*(-Mu*PLast-A*PLast+newEntries);
P_cur=PLast+.5*dt*(-Mu*PLast-A*PLast -Mu*PInt - A*PInt +2*newEntries );

%update solution
sols=[sols P_cur];

%two past solutions needed for BDF2
PLast2=PLast;
PLast=P_cur;
%BDF2
for i=1:nTimeSteps

    curTime=curTime+dt;
    timeOutput=[timeOutput;curTime];

    %update the 
    if(curTime-floor(curTime)==0 && curTime<=2023)    
        mm=curTime
        curYearData=PWHNewDiagData(PWHNewDiagData.Year==curTime,:);  
        [diagMu,diagSigma,newDiags,~]=...
            age_Distribution_MaxLikelihood(ages,[0;curYearData.Cases]);
            newEntries=sum(curYearData.Cases)*newDiags(ageMesh);
    end

    rhs=newEntries+2*M*PLast -.5*M*PLast2 - Mu*PLast;
    rhs(end)=0;
    P_cur=(1.5*M+A)\(rhs);    
    sols=[sols P_cur];    
    PLast2=PLast;
    PLast=P_cur;

end

plotTimeSeries(sols,timeOutput,ageMesh);
