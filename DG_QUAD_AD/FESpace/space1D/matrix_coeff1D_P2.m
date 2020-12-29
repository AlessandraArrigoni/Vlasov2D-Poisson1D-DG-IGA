%--------------------------------------------------------------------
% PURPOSE:
% 
% This routine computes the coefficients of the P2 Lagrange shape 
% functions 1D on the reference interval [0,1]
%
%--------------------------------------------------------------------

function [coef]= matrix_coeff1D_P2(nln)


%      x^2       x     c 
M=[
        0       0      1   % (0)
        1/4     1/2    1   % (0.5)
        1        1     1   % (1)
];

for i=1:nln
    f=zeros(nln,1);
    f(i)=1;
    coef(:,i)=M\f;
end