


function [dmdOutput,timeGrid] = dmd_Evaluation_Tool(apply_DMD,initialValue,initialTime,dt,nTimeSteps)
%Get an output by feeding the apply_DMD function here, together with the 
%initial value, initial time, time step, and number of time steps.

	dmdOutput=[initialValue];
    curTime=initialTime;
    timeGrid=[initialTime];
    
    
    for i=1:nTimeSteps
        curTime=curTime+dt;
        timeGrid=[timeGrid;curTime];
        currentOutput=apply_DMD(initialValue,(curTime-initialTime));
        dmdOutput=[dmdOutput currentOutput];
    end

end

