function [coef]= matrix_coeff_P2(nln)
%% [coef]= matrix_coeff_P2(nln)
% PURPOSE:
%
% This routine <omputes the coefficients of the P2 Lagrange shape 
% functions.
%     
% INPUT:
%       nln : (int) local dofs
% OUTPUT:
%       corff : (array real) coefficient of the linear expansion
%--------------------------------------------------------------------

f=zeros(nln,1);
coef_matrix=zeros(nln,nln);

%       x^2     y^2      xy       x     y      c 
M=[
        0       0        0       0      0       1      % (0,0)
       1/4      0        0      1/2     0       1      % (0.5,0)
        1       0        0       1      0       1      % (1,0)
       1/4     1/4      1/4     1/2    1/2      1      % (0.5,0.5)
        0       1        0       0      1       1      % (0,1)
        0      1/4       0       0     1/2      1      % (0,0.5)
];

for i=1:nln
    f=zeros(nln,1);
    f(i)=1;
    coef(:,i)=M\f;
end