function [empiricPDF] = age_Distribution_Empirical_FitDist(inputAge_Range,inputAge_Observations)
tic
%This takes an age-binned function and returns a distribution based on it.
%Matlab can't make distributions directly from frequency data, so we use
%the data from our CDC data to essentially "hack" this by making a fake
%sample set which follows the CDC data.
    
    randSamples=generateSamples(inputAge_Range,inputAge_Observations);
    empiricDistribution=fitdist(randSamples,'Kernel','Width',4);
    empiricPDF=@(x) pdf(empiricDistribution,x);

end

