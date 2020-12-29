function [coef]= matrix_coeff_P3(nln)
%% [coef]= matrix_coeff_P3(nln)
% PURPOSE:
%
% This routine <omputes the coefficients of the P3 Lagrange shape 
% functions.
%     
% INPUT:
%       nln : (int) local dofs
% OUTPUT:
%       corff : (array real) coefficient of the linear expansion
%--------------------------------------------------------------------

f=zeros(nln,1);
coef_matrix=zeros(nln,nln);

%     x^3  y^3  x^2y  y^2x   x^2   y^2      xy       x     y      c 
M=[
       0    0     0    0     0       0        0       0      0       1      % (0,0)
    1/64    0     0    0     1/16    0        0     1/4      0       1      % (1/4,0)
   27/64    0     0    0     9/16    0        0     3/4      0       1      % (3/4,0)
       1    0     0    0     1       0        0       1      0       1      % (1,0)
   27/64 1/64  9/64 3/64     9/16 1/16     3/16     3/4    1/4       1      % (3/4,1/4)
   1/64  27/64 3/64 9/64     1/16 9/16     3/16     1/4    3/4       1      % (1/4,3/4)
      0     1     0    0     0       1        0       0      1       1      % (0,1)
      0 27/64     0    0     0    9/16        0       0    3/4       1      % (0,3/4)
      0  1/64     0    0     0    1/16        0       0    1/4       1      % (0,1/4)
   1/27  1/27  1/27 1/27   1/9     1/9      1/9     1/3    1/3       1      % (1/3,1/3)
      
];


for i=1:nln
    f=zeros(nln,1);
    f(i)=1;
    coef(:,i)=M\f;
end