addpath('../Functions/');
load('../Data/PWHPrevData.mat');
load('../Data/PWHNewDiagData.mat');
load('../Data/lifeTable.mat');




%read data
firstYearData=PWHPrevData(PWHPrevData.Year==2008,:);
ages=[0;firstYearData.AgeGroup;100];
annPrev=[firstYearData.Cases];

covarFactor=.1;

%set some parameter values
minAge=1;
maxAge=101;


%this is the max likelihood fitting to the 2008 data to initialize the
%population
%[maxMu,maxSigma,popDist,~]=...
%    age_Distribution_MaxLikelihood(ages,[0;firstYearData.Cases]);
popDist=age_Distribution_Empirical_FitDist(ages,[0;firstYearData.Cases]);

%time step
dt=.25;

%age mesh initialization
ageMesh=linspace(minAge,maxAge,(maxAge-minAge)/dt)';

%lets us quickly (if roughly) evaluate the population sizes in each bin
[intMat,indSet]=buildAgeFunctionalMatrix(ages,ageMesh);

%system matrices
A=assemble_Age_Matrix(minAge,maxAge,dt);
M=assemble_Mass_Matrix(size(A,1),dt);
Mu=assemble_Mortality_Matrix(lifeTable,ageMesh);

%initial distribution
P_0=sum(firstYearData.Cases)*popDist(ageMesh);

%solution storage matrix
sols=[P_0];

%generate ensembles
ensemble=mvnrnd(P_0,.001*(P_0*P_0'),100)';
Q=ones(length(P_0),length(P_0));

%number of time steps
nTimeSteps=25/dt;


%beginningYear
yearStart=2008;
timeOutput=[yearStart];

%initialize RHS for first time step
curTime=yearStart+dt;
timeOutput=[timeOutput;curTime];
curYearData=PWHNewDiagData(PWHNewDiagData.Year==yearStart,:);  
%[diagMu,diagSigma,newDiags,~]=...
%    age_Distribution_MaxLikelihood(ages,[0;curYearData.Cases]);
newDiags=age_Distribution_Empirical_FitDist(ages,[0;curYearData.Cases]);
newEntries=sum(curYearData.Cases)*newDiags(ageMesh);
%"last" time step
PLast=ensemble;
newDiagData=[curYearData.Cases];
%needed for Heun's
PInt=zeros(size(ensemble));
P_cur=zeros(size(ensemble));

%first time step performed with with Heun's method
for j=1:size(ensemble,2)
    PInt(:,j)=PLast(:,j)+dt*(-Mu*PLast(:,j)-A*PLast(:,j)+newEntries);
    P_cur(:,j)=PLast(:,j)+.5*dt*(-Mu*PLast(:,j)-A*PLast(:,j) -Mu*PInt(:,j) - A*PInt(:,j) +2*newEntries );
end

%update solution
solsPreCor=[sols mean(P_cur,2)];
stateMean=mean(P_cur,2);

%EnKF (only active if dt=1)
if(curTime-floor(curTime)==0 )
    curYearData=PWHPrevData(PWHPrevData.Year==curTime,:);        
    ensembleFunctional=intMat*P_cur;        
    meanEnsemble=mean(ensembleFunctional,2);
    annPrev=[annPrev curYearData.Cases];
    fac1=ensembleFunctional-meanEnsemble*ones(100,1)';
    fac2=P_cur-stateMean*ones(100,1)';

    Pzz=(1/(size(ensemble,2)-1))*fac1*fac1'+diag(covarFactor*ones(7,1));
    Pxz=(1/(size(ensemble,2)-1))*fac2*fac1';

    K=Pxz*(Pzz\eye(7,7));
    
    P_cur=P_cur+K*([0;curYearData.Cases]*ones(100,1)'+mvnrnd(zeros(7,1),diag(covarFactor*ones(7,1)),100)' -ensembleFunctional);

    sols=[sols mean(P_cur,2)];

end

%two past solutions needed for BDF2
PLast2=PLast;
PLast=P_cur;

newDiagsByYear=[newEntries];
diagPDFsByYear=[newDiags(ageMesh)];

%begin BDF2 loop
for i=2:nTimeSteps

    curTime=curTime+dt;
    timeOutput=[timeOutput;curTime];

    %update the 
    if(curTime-floor(curTime)==0 && curTime<=2022)    
        mm=curTime
        curYearData=PWHNewDiagData(PWHNewDiagData.Year==curTime,:);  
        curYearPrevData=PWHPrevData(PWHPrevData.Year==curTime,:);  
%        [diagMu,diagSigma,newDiags,~]=...
%            age_Distribution_MaxLikelihood(ages,[0;curYearData.Cases]);
            newDiags=age_Distribution_Empirical_FitDist(ages,[0;curYearData.Cases]);
            newEntries=sum(curYearData.Cases)*newDiags(ageMesh);
            newDiagsByYear=[newDiagsByYear newEntries];
            diagPDFsByYear=[diagPDFsByYear newDiags(ageMesh)];
            annPrev=[annPrev curYearPrevData.Cases];
            newDiagData=[newDiagData curYearData.Cases];
            %nd=size(newDiagData)
            %anp=size(annPrev)
    elseif(curTime==2023)
        mm=curTime
         [dmdOp1,dmdOp2,apply_DMD]=dmd_Diag_Probability(newDiagsByYear,5,15,1);
         newEntries=apply_DMD(newDiagsByYear(:,1),(curTime-2008));
    elseif(curTime>2023)
        mm=curTime
        newEntries=apply_DMD(newDiagsByYear(:,1),(curTime-2008));
    end

    for j=1:size(ensemble,2)
        rhs=newEntries+2*M*PLast(:,j) -.5*M*PLast2(:,j) - Mu*PLast(:,j);
        rhs(end)=0;
        P_cur(:,j)=(1.5*M+A)\(rhs);    
        %sols=[sols P_cur];    
        %PLast2=PLast;
        %PLast=P_cur;
    end
    
    stateMean=mean(P_cur,2);
    solsPreCor=[solsPreCor stateMean];

    
    if(curTime-floor(curTime)==0 && curTime<=2022)    
        
        curYearData=PWHPrevData(PWHPrevData.Year==curTime,:);       
        ensembleFunctional=intMat*P_cur;        
        meanEnsemble=mean(ensembleFunctional,2); 

        fac1=ensembleFunctional-meanEnsemble*ones(100,1)';
        fac2=P_cur-stateMean*ones(100,1)';

        Pzz=(1/(size(ensemble,2)-1))*fac1*fac1'+diag(covarFactor*ones(7,1));
        Pxz=(1/(size(ensemble,2)-1))*fac2*fac1';
        
        K=Pxz*(Pzz\eye(7,7));

        P_cur=P_cur+K*([0;curYearData.Cases]*ones(100,1)' +mvnrnd(zeros(7,1),diag(covarFactor*ones(7,1)),100)' -ensembleFunctional);

    end
   
    PLast2=PLast;
    PLast=P_cur;
    
    sols=[sols mean(P_cur,2)];
    
end

%plotting function
plotTimeSeries(sols,timeOutput,ageMesh);
