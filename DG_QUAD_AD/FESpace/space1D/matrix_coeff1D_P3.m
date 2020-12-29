%--------------------------------------------------------------------
% PURPOSE:
% 
% This routine computes the coefficients of the P3 Lagrange shape 
% functions 1D on the reference interval [0,1]
%
%--------------------------------------------------------------------

function [coef]= matrix_coeff1D_P3(nln)


%       x^3      x^2      x     c 

M=[
        0        0       0      1  % (0)
        1/64     1/16    1/4    1  % (0.25)
        27/64    9/16    3/4    1  % (0.75)
        1        1        1     1  % ( 1)
];

for i=1:nln
    f=zeros(nln,1);
    f(i)=1;
    coef(:,i)=M\f;
end