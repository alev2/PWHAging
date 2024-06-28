function [muMax,sigmaMax,maxLogNormPDF, likelihoodFun] = age_Distribution_MaxLikelihood(inputAge_Range,inputAge_Obseravations)

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
logNormPDF= @(x,mu,sigma)...
    (1./(x*sigma*sqrt(2*pi))).*exp(-(log(x)-mu).^2./(2*sigma^2));

%the parameter space for max likelihood
xx=linspace(log(5),log(80),50);
mumu=linspace(log(1.05),log(2),50);

%this is an easy problem, so we get the max likelihood with a brute force
%approach. Bite me.
likelihoodFun=[];
for i=1:length(xx)
    for j=1:length(mumu)
        likelihood=0;
        for k=1:(length(inputAge_Range)-1)                
            likelihood=...
               likelihood...
                +log((integral(@(x)logNormPDF(x,xx(i),mumu(j)),inputAge_Range(k),(inputAge_Range(k+1)-1))^ageObservations(k)));
        end
        likelihoodFun=[likelihoodFun; exp(xx(i)) exp(mumu(j)) likelihood];
    end
end

%the max likeliood
[~,maxInd]=max(likelihoodFun(:,3));

%return the parameters and the distributions
muMax=likelihoodFun(maxInd,1);
sigmaMax=likelihoodFun(maxInd,2);

maxLogNormPDF= @(x)...
    (1./(x*log(sigmaMax)*sqrt(2*pi))).*exp(-(log(x)-log(muMax)).^2./(2*log(sigmaMax)^2));

end

