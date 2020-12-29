%--------------------------------------------------------------------
% PURPOSE:
% 
% This routine computes the coefficients of the P1 Lagrange shape 
% functions.
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [coef]= matrix_coeff_P1(nln);

f=zeros(nln,1);
coef_matrix=zeros(nln,nln);


%       x        y      c 
M=[
        0       0       1      % (0,0)
        1       0       1      % (1,0)
        0       1       1      % (0,1)
];

for i=1:nln
    f=zeros(nln,1);
    f(i)=1;
    coef(:,i)=M\f;
end