addpath('../Functions/');
load('../Data/PWHPrevData.mat');
load('../Data/PWHNewDiagData.mat');
load('../Data/PWHDeathData.mat');
load('../Data/lifeTable.mat');

%this is the code that seeks to define the age-dependent PWH mortality in
%time (Iacomini & Sgattoni project)


%read data
firstYearData=PWHPrevData(PWHPrevData.Year==2008,:);
ages=[0;firstYearData.AgeGroup;100];
%agesK=[

annPrev=[firstYearData.Cases];
nParticles=50;
covarFactor=1e0;%1e-7;
ensembleFactor=1e-2;
%set some parameter values
minAge=1;
maxAge=101;

numFilterPasses=2;


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
%Mu=assemble_Mortality_Matrix(lifeTable,ageMesh);

%initial distribution
P_0=sum(firstYearData.Cases)*popDist(ageMesh);

%solution storage matrix
solsState=[P_0];
%generate ensembles
muBase=diag(assemble_Mortality_Matrix(lifeTable,ageMesh));
sols=[muBase];
ensemble=mvnrnd(muBase,ensembleFactor*(muBase*muBase'),nParticles)';
Q=ones(length(muBase),length(muBase));

%number of time steps
nTimeSteps=15/dt;


%beginningYear
yearStart=2008;
timeOutput=[yearStart];

%initialize RHS for first time step
mm=yearStart;
curTime=yearStart+dt;
timeOutput=[timeOutput;curTime];
curYearData=PWHNewDiagData(PWHNewDiagData.Year==yearStart,:);  
%[diagMu,diagSigma,newDiags,~]=...
%    age_Distribution_MaxLikelihood(ages,[0;curYearData.Cases]);
newDiags=age_Distribution_Empirical_FitDist(ages,[0;curYearData.Cases]);
newEntries=sum(curYearData.Cases)*newDiags(ageMesh);
%"last" time step
PLast=P_0*ones(size(ensemble,2),1)';
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
    Mu=diag(ensemble(:,j));
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
    
    Pzz=(1/(size(ensemble,2)-1))*fac1*fac1'+diag(covarFactor*ones(7,1));
    Pxz=(1/(size(ensemble,2)-1))*fac2*fac1';

    K=Pxz*(Pzz\eye(7,7));
    
    ensemble=ensemble+K*([0;curYearData.Cases]*ones(nParticles,1)'+mvnrnd(zeros(7,1),diag(covarFactor*ones(7,1)),nParticles)' -ensembleFunctional);

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
while curTime<=(yearStart+nTimeSteps*dt)
    
    curTime=curTime+dt
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
            deceasedCur=zeros(size(PLast));
         PLast2=curYrLast2;
         PLast=curYrLast;
    elseif(curTime==2023)
        mm=curTime
         [dmdOp1,dmdOp2,apply_DMD]=dmd_Diag_Probability(newDiagsByYear,5,15,1);
         newEntries=apply_DMD(newDiagsByYear(:,1),(curTime-2008));
         deceasedCur=zeros(size(PLast));
    elseif(curTime>2023)
        mm=curTime
        newEntries=apply_DMD(newDiagsByYear(:,1),(curTime-2008));
        deceasedCur=zeros(size(PLast));
    end

    for j=1:size(ensemble,2)
        %%%REEXAMINE WHETHER IT'S ACTUALLY THE BEST IDEA TO ADVANCE THE
        %%%MORTALITY EXPLICITLY... SOMEHOW I DOUBT IT.
        
        %Mu=assemble_Mortality_Matrix([lifeTable(:,1) ensemble(:,j)],ageMesh);
        %Mu=real(Mu);
        Mu=diag(ensemble(:,j));
        
        rhs=newEntries+2*M*PLast(:,j) -.5*M*PLast2(:,j);% - Mu*PLast(:,j);
        rhs(end)=0;
        
        P_cur(:,j)=(1.5*M+A+Mu)\(rhs);    
        
        deceasedCur(:,j)= deceasedCur(:,j) + dt*Mu*P_cur(:,j);
        %deceasedCur(:,j)= deceasedCur(:,j) + dt*Mu*PLast(:,j);        

        %sols=[sols P_cur];    
        %PLast2=PLast;
        %PLast=P_cur;
    end
    stateMean=mean(P_cur,2);
    solsStatePreCor=[solsStatePreCor stateMean];
    
    PLast2=PLast;
    PLast=P_cur;
   
    if((curTime+dt)-floor(curTime+dt)==0 && (curTime+dt)<=2022)    

        filterPass=filterPass+1;
        curState=[curTime filterPass]
        mortRateMean=mean(ensemble,2);

        curYearData=PWHDeathData(PWHDeathData.Year==(curTime+dt),:)          
        ensembleFunctional=intMat*deceasedCur;        
        meanEnsemble=mean(ensembleFunctional,2);
    
        annPrev=[annPrev curYearData.Cases];
        fac1=ensembleFunctional-meanEnsemble*ones(nParticles,1)';
        fac2=ensemble-mortRateMean*ones(nParticles,1)';
    
        Pzz=(1/(size(ensemble,2)-1))*fac1*fac1'+diag(covarFactor*ones(7,1));
        Pxz=(1/(size(ensemble,2)-1))*fac2*fac1';

        K=Pxz*(Pzz\eye(7,7));
    
        ensemble=ensemble+K*([0;curYearData.Cases]*ones(nParticles,1)'+mvnrnd(zeros(7,1),diag(covarFactor*ones(7,1)),nParticles)' -ensembleFunctional);

    
        %deceasedCur=zeros(size(ensemble));
        if(filterPass<=numFilterPasses)
           curTime=mm;
           deceasedCur=zeros(size(ensemble));
           PLast2=curYrLast2;
           PLast=curYrLast;
        else            
           deceasedCur=zeros(size(ensemble));
           filterPass=0;
           solsState=[solsState mean(P_cur,2)];
           solTimes=[solTimes; curTime+dt];
           sols=[sols mean(ensemble,2)];
           deathsByYear=[deathsByYear meanEnsemble];
           deathsByYear_Correction=[deathsByYear_Correction intMat*(diag(mean(ensemble,2))*stateMean)];
           deathsByYear_Real=[deathsByYear_Real curYearData.Cases];
           curYrLast2=PLast;
           curYrLast=P_cur;
        end
        
    end
   
    
    
end

%plotting function
%plotTimeSeries(solsState,timeOutput,ageMesh);
generateInvKF_BarPlot(deathsByYear_Real,deathsByYear(2:end,:),solTimes(2:end),ageMesh,sols);