function [intMat,indSet] = buildAgeFunctionalMatrix(ages,ageMesh)
%this builds aa matrix to evaluate the population sizes inside each
%population bin. The age mesh is assumed uniform.

    da=ageMesh(2)-ageMesh(1);
    
    indSet=1;
    intMat=zeros(size(ages,1)-1,size(ageMesh,1));
    
    for i=2:length(ages)
       
        
       [~,b]=min(abs(ageMesh-ages(i)));       
       indSet=[indSet;b]; 
        
       intMat(i-1,indSet(i-1):indSet(i))=da;
       
    end
%    outputArg2 = inputArg2;


end

