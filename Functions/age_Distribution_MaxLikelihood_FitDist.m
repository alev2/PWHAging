function [muMax,sigmaMax,maxLogNormPDF, likelihoodFun] = age_Distribution_MaxLikelihood_FitDist(inputAge_Range,inputAge_Obseravations)
tic
%This function takes aan n-dimensional vector formatted like this:
%  [0
%   13
%   24
%   35
%   45
%   55
%   100]
% means that we group into bins of 0-13, 13-24, 25-34, and so forth
%
%   Corresponding probabilities vector:
%   [0
%    0.05
%   ...
% means probability of 0 in the first group, .05 in the second, so forth.


%get the "PDF" of the age distribution
ageObservations=round(100*inputAge_Obseravations/sum(inputAge_Obseravations));

%function to be evaluated
%logNormPDF= @(x,mu,sigma)...
%    (1./(x*sigma*sqrt(2*pi))).*exp(-(log(x)-mu).^2./(2*sigma^2));
logNormPDF= @(x,mu,sigma)lognpdf(x,mu,sigma);
%gammaPDF= @(x,mu,sigma)gampdf(x,mu,sigma);

%the parameter space for max likelihood
xx=linspace(log(5),log(80),150);
mumu=linspace(log(1.05),log(2),150);

%xx=linspace(5,80,150);
%mumu=linspace(1.0,10,150);


%this is an easy problem, so we get the max likelihood with a brute force
%approach.
likelihoodFun=[];
for i=1:length(xx)
    for j=1:length(mumu)
        likelihood=evaluateLikelihood(xx(i),mumu(j),inputAge_Range,ageObservations);
        likelihoodFun=[likelihoodFun; exp(xx(i)) exp(mumu(j)) likelihood];
%        likelihoodFun=[likelihoodFun; xx(i) mumu(j) likelihood];
    end
end






%the max likeliood
[~,maxInd]=max(likelihoodFun(:,3));

%return the parameters and the distributions
muMax=likelihoodFun(maxInd,1);
sigmaMax=likelihoodFun(maxInd,2);

   
function likeEval=evaluateLikelihood(input_x,input_mu,input_Age,ageObs)
    likeEval=0;
        for kk=1:(length(input_Age)-1)                
            likeEval=...
                likeEval+...
                ...ageObs(kk)*log( (gamcdf(input_Age(kk+1),input_x,input_mu)-gamcdf(input_Age(kk),input_x,input_mu)));
                ageObs(kk)*log( (logncdf(input_Age(kk+1),input_x,input_mu)-logncdf(input_Age(kk),input_x,input_mu)));
       end    
end
maxLogNormPDF= @(x)lognpdf(x,log(muMax),log(sigmaMax));
%maxLogNormPDF= @(x)gampdf(x,(muMax),(sigmaMax));

toc        
end

