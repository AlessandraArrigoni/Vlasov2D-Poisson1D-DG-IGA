%--------------------------------------------------------------------
% PURPOSE:
% 
% This routine computes the coefficients of the Q1 Lagrange shape 
% functions on the reference element [-1,1]x[-1,1]
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [coef]= matrix_coeff_Q1(nln)

f=zeros(nln,1);
coef_matrix=zeros(nln,nln);


%       x        y      xy    c 
M=[
       -1        1      -1     1      % (-1, 1)
       -1       -1       1     1      % (-1,-1)
        1       -1      -1     1      % ( 1,-1)
        1        1       1     1      % ( 1, 1)
];

for i=1:nln
    f=zeros(nln,1);
    f(i)=1;
    coef(:,i)=M\f;
end