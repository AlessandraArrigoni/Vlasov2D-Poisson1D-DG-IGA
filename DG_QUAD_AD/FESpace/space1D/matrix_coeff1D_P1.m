%--------------------------------------------------------------------
% PURPOSE:
% 
% This routine computes the coefficients of the P1 Lagrange shape 
% functions 1D on the reference interval [0,1]
%
%--------------------------------------------------------------------

function [coef]= matrix_coeff1D_P1(nln)


%       x        c 
M=[
        0        1    % (0)
        1        1    % (1)
];

for i=1:nln
    f=zeros(nln,1);
    f(i)=1;
    coef(:,i)=M\f;
end