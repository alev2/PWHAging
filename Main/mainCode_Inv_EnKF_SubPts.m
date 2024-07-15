addpath('../Functions/');
load('../Data/PWHPrevData_FinerStrat.mat');
load('../Data/PWHNewDiagData_FinerStrat.mat');
load('../Data/PWHDeathData_FinerStrat.mat');
%load('../Data/PWHPrevData.mat');
%load('../Data/PWHNewDiagData.mat');
%load('../Data/PWHDeathData.mat');

load('../Data/lifeTable.mat');

%this is the code that seeks to define the age-dependent PWH mortality in
%time (Iacomini & Sgattoni project).

%this seeks to solve the inverse kalman filter problem.


%read data
firstYearData=PWHPrevData(PWHPrevData.Year==2008,:);
firstYearDeathData=PWHDeathData(PWHDeathData.Year==2008,:);
ages=[0;firstYearData.AgeGroup;100];

interpType='pchip';
%helpful for printing

%good
%covarfactor 1e2, ensfactor 1e-3, tikh 1e-5, nParticles 50, filtPasses100

%

annPrev=[firstYearDeathData.Cases];
nParticles=50;
covarFactor=1e5;%1e-7;
ensembleFactor=.05;
filterWidth=1;
skipStep=60;
tikh=1e3;
lowFactor=0;%1e-6;%1e-8;%1e2;

%time step
dt=.25;


%set some parqameter values
minAge=1;
maxAge=101;
ageMesh=linspace(minAge,maxAge,1+((maxAge-minAge)/dt))';

numFilterPasses=25;


%this is the max likelihood fitting to the 2008 data to initialize the
%population
%[maxMu,maxSigma,popDist,~]=...
%    age_Distribution_MaxLikelihood(ages,[0;firstYearData.Cases]);

P_0=[];

for kk=1:nParticles
    
    popDist=age_Distribution_Empirical_FitDist(ages,[0;firstYearData.Cases]);
    P_0=[P_0 sum(firstYearData.Cases)*popDist(ageMesh)];
    
end

ages=[0;firstYearDeathData.AgeGroup;100];
nBins=size(ages,1)-1;


%age mesh initialization
rrr=0;
%ageMesh2=[ageMesh(1:skipStep:end)];
%ageMesh2=[ages;101];
%rrr=1;
%if (ageMesh2(end)~=101)
%   ageMesh2=[ageMesh2;101]; 
%   rrr=1;
%end

agesKK= [ages(1);8;ages(2:end-1); ages(end)];
ageMesh2=[];
ageMesh2Inds=[];
for kk=1:length(agesKK)
   [~,bb]=min(abs(ageMesh-agesKK(kk)) );
   ageMesh2=[ageMesh2;ageMesh(bb)];
   ageMesh2Inds=[ageMesh2Inds;bb];
end

ageMesh2=[ageMesh2;101];
ageMesh2Inds=[ageMesh2Inds;size(ageMesh,1)];
%ageMesh2=[ageMesh2(1);.5*(ageMesh2(1)+ageMesh2(2));ageMesh2(2:end)];
%lets us quickly (if roughly) evaluate the population sizes in each bin
[intMat,indSet]=buildAgeFunctionalMatrix(ages,ageMesh);

initialDeathVariability=round(intMat*(diag(interpn(ageMesh2,muBase,ageMesh,interpType))*P_0))

%system matrices
A=assemble_Age_Matrix(minAge,maxAge,dt);
M=assemble_Mass_Matrix(size(A,1),dt);
%Mu=assemble_Mortality_Matrix(lifeTable,ageMesh);

%initial distribution


penaltySize=22-1;
rankR=100;

zz1=zeros(size(intMat,1),penaltySize);
zz2=zeros(penaltySize,size(intMat,2));
funcMat=[intMat zz1;  zz2 eye(penaltySize,penaltySize)];

%solution storage matrix
solsState=[P_0];
%generate ensembles
%lifeTable(end,2)=lifeTable(end,2)-eps;
muBaseKK=diag(assemble_Mortality_Matrix(lifeTable,ageMesh));
muBase=muBaseKK(ageMesh2Inds);

%if (rrr==1)
%   muBase=[muBase;muBaseKK(end)]; 
%end


sols=[muBase];
Q=ones(length(muBase),length(muBase));


%XX=makeLocalCovarianceMatrix(length(muBase));
%XX=diag((muBase))*XX*diag((muBase));
%ensemble=(mvnrnd(muBase,ensembleFactor*XX,nParticles)');

ensemble=mvnrnd(muBase,ensembleFactor*diag(muBase.*muBase),nParticles)';
%ensemble=mvnrnd(muBase,ensembleFactor*diag(ones(size(muBase,1),1)),nParticles)';
%ensemble=mvnrnd(muBase,ensembleFactor*diag(muBase.*muBase),nParticles)';
%ensemble=ensemble+.1*randn(size(ensemble)).*ensemble;
%ensemble=muBase*ones(nParticles,1)'+((.025)*muBase*rand(nParticles,1)')
%ensemble=log(MvLogNRand(muBase,ensembleFactor*(muBase),nParticles,diag(muBase)*ensembleFactor)');

ensMin=min(min(ensemble,[],2))


%number of time steps
nTimeSteps=1+14/dt;


%beginningYear
yearStart=2009;
timeOutput=[yearStart];

%initialize RHS for first time step
filterTime=yearStart;
curTime=yearStart+dt;
timeOutput=[timeOutput;curTime];
curYearData=PWHNewDiagData(PWHNewDiagData.Year==yearStart,:);  
%[diagMu,diagSigma,newDiags,~]=...
%    age_Distribution_MaxLikelihood(ages,[0;curYearData.Cases]);
newDiags=age_Distribution_Empirical_FitDist(ages,[0;curYearData.Cases]);
newEntries=sum(curYearData.Cases)*newDiags(ageMesh);
%"last" time step
PLast=P_0;%*ones(size(ensemble,2),1)';
newDiagData=[curYearData.Cases];
%needed for Heun's
PInt=zeros(size(PLast));
P_cur=zeros(size(PLast));

%keeps track of current year deceased
deceasedCur=zeros(size(PLast));
deathsByYear=[];
deathsByYear_Correction=[];
deathsByYear_Real=[];
%first time step performed with with Heun's method
for j=1:size(ensemble,2)
    %Mu=assemble_Mortality_Matrix([lifeTable(:,1) ensemble(:,j)],ageMesh);
	Mu=diag(interpn(ageMesh2,ensemble(:,j),ageMesh,interpType));
    
    PInt(:,j)=PLast(:,j)+dt*(-Mu*PLast(:,j)-A*PLast(:,j)+newEntries);
    P_cur(:,j)=PLast(:,j)+.5*dt*(-Mu*PLast(:,j)-A*PLast(:,j) -Mu*PInt(:,j) - A*PInt(:,j) +2*newEntries );
    
    deceasedCur(:,j) = .5*dt*(Mu*(PLast(:,j)+PInt(:,j)));
    
end

%update solution
solsStatePreCor=[solsState mean(P_cur,2)];
decPreCor=[mean(deceasedCur,2)];
stateMean=mean(P_cur,2);
mortRateMean=mean(ensemble,2);
decMean=decPreCor;


%EnKF (only active if dt=1)
if(curTime-floor(curTime)==0 )
    
    curYearData=PWHDeathData(PWHDeathData.Year==curTime,:);        
    ensembleFunctional=intMat*deceasedCur;        

    meanEnsemble=mean(ensembleFunctional,2);
    annPrev=[annPrev curYearData.Cases];
    fac1=ensembleFunctional-meanEnsemble*ones(nParticles,1)';
    fac2=ensemble-mortRateMean*ones(nParticles,1)';
    
    Pzz=(1/(size(ensemble,2)-1))*(fac1*fac1')+diag(covarFactor*[ones(nBins,1)]);
%    Pzz=(1/(size(ensemble,2)-1))*(fac1*fac1')+(covarFactor*makeLocalCovarianceMatrix(nBins));
    Pxz=(1/(size(ensemble,2)-1))*fac2*fac1';

    K=Pxz*(Pzz\eye(nBins,nBins));
   
    %ensemble=ensemble+K*( max( [1;curYearData.Cases]*ones(nParticles,1)'+ mvnrnd(zeros(nBins,1),(covarFactor*diag([ ones(nBins,1)])),nParticles)',1) -ensembleFunctional);
    ensemble=ensemble+K*( [1;curYearData.Cases]*ones(nParticles,1)'+ mvnrnd(zeros(nBins,1),(covarFactor*diag([ ones(nBins,1)])),nParticles)' -ensembleFunctional);
    %ensemble=ensemble+K*( max( [0;curYearData.Cases]*ones(nParticles,1)'+mvnrnd(zeros(nBins,1),(covarFactor*diag([ones(nBins,1)])),nParticles)',0) -ensembleFunctional);
    %ensemble=ensemble+K*([0;curYearData.Cases]*ones(nParticles,1)'+mvnrnd(zeros(nBins,1),(covarFactor*makeLocalCovarianceMatrix(nBins)),nParticles)' -ensembleFunctional);
    %ensemble=movmean(ensemble+K*([0;curYearData.Cases]*ones(nParticles,1)' -ensembleFunctional),filterWidth);
    %ensemble=max(ensemble,muBase*ones(nParticles,1)');
    
    sols=[sols mean(ensemble,2)];
    deathsByYear=[deathsByYear meanEnsemble];
    deathsByYear_Real=[deathsByYear_Real curYearData.Cases];
    %deceasedCur=zeros(size(ensemble));

end

%two past solutions needed for BDF2
PLast2=PLast;
PLast=P_cur;

curYrLast2=PLast;
curYrLast=P_cur;

newDiagsByYear=[newEntries];
diagPDFsByYear=[newDiags(ageMesh)];

filterPass=0;
%begin BDF2 loop
solTimes=yearStart;
ensembleLast=ensemble;
convergenceHistories=[];

while curTime<=(yearStart+nTimeSteps*dt)
    
    curTime=curTime+dt;
    timeOutput=[timeOutput;curTime];

    %update the data at the correct years.
    if(curTime-floor(curTime)==0 && curTime<2023)    
        %keep this so we can reset time counter when we do filter passes
        filterTime=curTime;

        %current year diagnosis and prevalence data.
        curYearData=PWHNewDiagData(PWHNewDiagData.Year==curTime,:);  
        curYearPrevData=PWHPrevData(PWHPrevData.Year==(curTime-1),:);  

        %old method - not useful and not used anymore.
        %       [diagMu,diagSigma,newDiags,~]=...
        %            age_Distribution_MaxLikelihood(ages,[0;curYearData.Cases]);
        
        %current year new diagnoses and entries
        newDiags=age_Distribution_Empirical_FitDist(ages,[0;curYearData.Cases]);
        newEntries=sum(curYearData.Cases)*newDiags(ageMesh);

        %initialize certain current step variables.
        deceasedCur=zeros(size(PLast));
        %at each filter pass we need to reset these.
        PLast2=curYrLast2;
        PLast=curYrLast;        
        
        if (filterPass==0)
            %misc outputs for post-processing
            newDiagsByYear=[newDiagsByYear newEntries];
            diagPDFsByYear=[diagPDFsByYear newDiags(ageMesh)];
            newDiagData=[newDiagData curYearData.Cases];
        end
        
    elseif(curTime==2023)
        filterTime=curTime;
         [dmdOp1,dmdOp2,apply_DMD]=dmd_Diag_Probability(newDiagsByYear,5,size(newDiagsByYear,2),1);
         newEntries=apply_DMD(newDiagsByYear(:,1),(curTime-2009));
         deceasedCur=zeros(size(PLast));
    elseif(curTime>2023)
        filterTime=curTime;
        newEntries=apply_DMD(newDiagsByYear(:,1),(curTime-2009));
        deceasedCur=zeros(size(PLast));
    end

    for j=1:size(ensemble,2)
        %%This solves the solution state for each ensemble member.
        
        %build mortality matrix
        Mu=diag(interpn(ageMesh2,ensemble(:,j),ageMesh,interpType));

        %assemble RHS and apply boundary conditions
        rhs=newEntries+2*M*PLast(:,j) -.5*M*PLast2(:,j);
        rhs(end)=0;

        %solution state for each ensemble member.
        P_cur(:,j)=(1.5*M+A+Mu)\(rhs);    
        
        %keep track of current deceased state for each ensemble member
        %(pre-binning)
        deceasedCur(:,j)= deceasedCur(:,j) + dt*Mu*P_cur(:,j);
    end
    
    %mean of solution state
    stateMean=mean(P_cur,2);
    
    %solution state BEFORE Kalman correction
    solsStatePreCor=[solsStatePreCor stateMean];
    
    %Update last steps for BDF2 advancement
    PLast2=PLast;
    PLast=P_cur;

    
    %EnKF procedure. Only do this once per year.
    if((curTime+dt)-floor(curTime+dt)==0 && (curTime+dt)<=2023)    

        %update filter pass
        filterPass=filterPass+1;
        
        %the current year and filter pass
        stateMessage=strcat('Year-end: ',num2str(floor(curTime),'%.0f'),'. Filter pass: ',num2str(filterPass),'.');
        disp(stateMessage);
        
        %mean of the Kalman ensemble for the mortality rate.
        mortRateMean=mean(ensemble,2);
        
        %this is useless at the moment.
        ensembleK=[];
        
        for kk=1:size(ensemble,2)
           curVec=interpn(ageMesh2,ensemble(:,kk),ageMesh,interpType);
           ensembleK=[ensembleK curVec];
        end
        
        %get the real data from the current simulation year.
        curYearData=PWHDeathData(PWHDeathData.Year==floor(curTime),:);          

        %evaluate ensemble functional.
        ensembleFunctional=intMat*[deceasedCur];
        ensembleFunctional=[ensembleFunctional];%; (diag([ones(390,1);zeros(10,1)])*(diff(ensembleK)/dt))];
        meanEnsemble=mean(ensembleFunctional,2);
            
        %annPrev=[annPrev curYearData.Cases];
        
        %build kalman matrices
        fac1=ensembleFunctional-meanEnsemble*ones(nParticles,1)';
        fac2=ensemble-mortRateMean*ones(nParticles,1)';
        
        %with tikhinov regularization on derivative of curve.
        Pzz=(1/(size(ensemble,2)-1))*(fac1*fac1')+diag([covarFactor*ones(nBins,1)]);%; tikh*ones(size(ensembleK,1)-1,1)]);
        
        %not necessary (for low-rank apx of Pzz).
        %[U,S,V]=svd(Pzz);    
        %UApx=U(:,1:rankR);
        %SApx=S(:,1:rankR);
        %SApx=diag(1./diag(SApx));
        %VApx=V(:,1:rankR);  
        %KFac=VApx*SApx*UApx';
        
        Pxz=(1/(size(ensemble,2)-1))*fac2*fac1';
        
        %K=Pxz*(Pzz\eye(nBins+size(ensembleK,1)-1,nBins+size(ensembleK,1)-1));
        K=Pxz*(Pzz\eye(nBins,nBins));
        %K=Pxz*KFac;
        
        %here we evalate the KF
        topBlock=[0;curYearData.Cases]*ones(nParticles,1)'+ mvnrnd(zeros(nBins,1),(covarFactor*diag([ ones(nBins,1)])),nParticles)';           
        %bottomBlock=mvnrnd( zeros(size(ensembleK,1)-1,1), lowFactor*diag(abs(diff(muBaseKK)/dt)) ,nParticles)';   
        MM=[topBlock];%;bottomBlock];
        ensemble=(ensemble+K*( MM - ensembleFunctional));
                  
        
        %In this section of the code, we update the filter.

        %CASE I: if we're under the number of filter passes, we do another
        %iteration of the current year.
        if(filterPass<=numFilterPasses)           
           %reset clock to current year
           curTime=filterTime;
           %reset deceased functional for current iteration
           deceasedCur=zeros(size(PLast));
           %last two state time steps reset to last steps from the prev yr
           PLast2=curYrLast2;
           PLast=curYrLast;           
           
           %track convergence
           if filterPass==1
              refMean=norm(mean(ensemble,2));
              annualConvergence=[];
           else
              convergenceCriterion=norm(mean(ensemble-ensembleLast,2))/refMean
              annualConvergence=[annualConvergence;convergenceCriterion];
           end
           
           ensembleLast=ensemble;
           
        %CASE II: if we have reached the number of filter passes, push
        %forward!
        else     
           %Reset deceased function and filter pass counter
           deceasedCur=zeros(size(PLast));
           filterPass=0;
           %save some output variables
           convergenceHistories=[convergenceHistories annualConvergence];
           solsState=[solsState mean(P_cur,2)];
           solTimes=[solTimes; curTime+dt];
           sols=[sols mean(ensemble,2)];
           deathsByYear=[deathsByYear meanEnsemble];
           ensNew=interpn(ageMesh2,mean(ensemble,2),ageMesh,interpType);
           deathsByYear_Correction=[deathsByYear_Correction intMat*(diag(ensNew)*stateMean)];
           deathsByYear_Real=[deathsByYear_Real curYearData.Cases];
           %Keep last two state values from this year, necessary for each filter
           %pass of the next year.
           curYrLast2=PLast;
           curYrLast=P_cur;
        end
        
    end
   
    
    
end


convergenceHistories=[solTimes(1:end-1)';convergenceHistories];

deathsByYear_Real=[deathsByYear_Real(1:4,:);sum(deathsByYear_Real(5:6,:),1);sum(deathsByYear_Real(7:end,:),1)];
deathsByYear=[deathsByYear(2:5,:);sum(deathsByYear(6:7,:),1);sum(deathsByYear(8:end,:),1)];


%plotting function
%plotTimeSeries(solsState,timeOutput,ageMesh);
generateInvKF_BarPlot_SubPts(deathsByYear_Real,deathsByYear,solTimes(1:(end-1)),ageMesh,sols,ageMesh2,interpType);