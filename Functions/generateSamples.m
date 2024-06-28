function [randSamples] = generateSamples(ageBins,observationCounts)

%This generates a set of fake observations in order to produce a sample
%that follows the age distribtions reported by CDC in the age data. This
%way we can directly fit a distribution from data. The alternatives are
%much less appealing.
    
    observationCountsNorm=round(500*observationCounts/sum(observationCounts));

 randSamples=[];

 for i=1:(length(ageBins)-1 )  
    
     curBin=ageBins(i)+(ageBins(i+1)-ageBins(i))*rand(observationCountsNorm(i),1);     
     randSamples=[randSamples;curBin];
 end        
    
    

end

