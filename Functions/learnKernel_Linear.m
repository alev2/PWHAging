function [K] = learnKernel_Linear(inputQuantity,outputQuantity,rank)
%This function provides a straightforward way of learning a state-dependent
%kernel as a linear function.
%Given two quantities (for example prevalence and annual diagnoses), it
%constructs the linear mapping that produces the QOI with a M-P inv

    [U,S,V]=svd(inputQuantity);
    
    UApx=U(:,1:rank);
    SApx=S(:,1:rank);
    SApx=diag(1./diag(SApx));
    VApx=V(:,1:rank);  
    
    K=outputQuantity*VApx*SApx*UApx';

end

