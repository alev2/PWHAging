function [dmdOperator,dmdOperator2,apply_DMD] = dmd_Diag_Probability(diagPDFsByYear, dmdRank, numTrainPoints,dt )

%This uses a DMD approach to obtain the time-extrapolated diagnosis age
%curve in time.

    X1=diagPDFsByYear(:,1:numTrainPoints-1);
    X2=diagPDFsByYear(:,2:numTrainPoints);

    [U1,S1,V1]=svd(X1);    
    [U2,S2,V2]=svd(X2);

    U1Apx=U1(:,1:dmdRank);
    S1Apx=S1(:,1:dmdRank);
    S1Apx=diag(1./diag(S1Apx));
    V1Apx=V1(:,1:dmdRank);

    dmdOperator=U2*S2*V2'*V1Apx*S1Apx*U1Apx';
    
    dmdOperator2=U1Apx'*X2*V1Apx*S1Apx;
    
    [W,Gamma]=eig(dmdOperator2);

    Psi=X2*V1Apx*S1Apx*W;

    
    apply_DMD=@(u,t)...
        real(Psi*exp((log(Gamma)*t)/dt)*pinv(Psi)*u);
end

