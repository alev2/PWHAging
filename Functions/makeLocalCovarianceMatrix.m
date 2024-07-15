function [X] = makeLocalCovarianceMatrix(N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    X=eye(N,N);
    
    for k=1:N
       
        X=X+diag(ones(N-k,1),k)*(2^(-k));
        X=X+diag(ones(N-k,1),-k)*(2^(-k));
        
    end

end

