function [A] = assemble_Age_Matrix(a,b,da)
%This returns the basic matrix for the aging matrid system.
%It corresponds to a three-term upwinded scheme.

    
%matrix size and such
    N=1+(b-a)/da;
    A=zeros(N,N);

%We can treat the inflow boundary naively here, because the population at
%the lowest ages is going to be zero pretty much the whole time. So
% 1.5u_1 \approx 1.5u_1  - 2u_{0} + .5u_{-1} 
%and so forth...

   A=A+...
       diag(ones(length(diag(A,0)),1)*3./(2*da),0)+...
       diag(ones(length(diag(A,-1)),1)*-4./(2*da),-1)+...
       diag(ones(length(diag(A,-2)),1)*1./(2*da),-2);

%     A=A+...
%         diag(ones(length(diag(A,0)),1)*1./(da),0)+...
%         diag(ones(length(diag(A,-1)),1)*-1./(da),-1);%+...
    
  
end

