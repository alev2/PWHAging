function [M] = assemble_Mass_Matrix(N,dt)
%This matrix is responsible for advancing      

    M=eye(N)./(1*dt);
   
    M(N,N)=0';
end

