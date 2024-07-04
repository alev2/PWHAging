function [randSamples] = generateSamples(ageBins,observationCounts)

%This generates a set of fake observations in order to produce a sample
%that follows the age distribtions reported by CDC in the age data. This
%way we can directly fit a distribution from data. The alternatives are
%much less appealing.
    
    observationCountsNorm=round(500*observationCounts/sum(observationCounts));

 randSamples=[];

 for i=1:(length(ageBins)-1 )  
    
    %if(i==(length(ageBins)-1))
      % curBin=ageBins(i)+(ageBins(i+1)-ageBins(i))*random('Normal',.1,.25,observationCountsNorm(i),1);     
    %  curBin=random('Poisson',ageBins(i)-8,observationCountsNorm(i),1);%,ageBins(i)+(ageBins(i+1)-ageBins(i))*random('Normal',.1,.25,observationCountsNorm(i),1);     
    %elseif(i==(length(ageBins)-2))
        %curBin=ageBins(i)+(ageBins(i+1)-ageBins(i))*random('Normal',.1,.25,observationCountsNorm(i),1);     
    %  curBin=random('Poisson',ageBins(i)+2,observationCountsNorm(i),1);%,ageBins(i)+(ageBins(i+1)-ageBins(i))*random('Normal',.1,.25,observationCountsNorm(i),1);     
    %else    
          curBin=ageBins(i)+(ageBins(i+1)-ageBins(i))*rand(observationCountsNorm(i),1);     
    %end
     randSamples=[randSamples;curBin];
 end        
    
    

end

