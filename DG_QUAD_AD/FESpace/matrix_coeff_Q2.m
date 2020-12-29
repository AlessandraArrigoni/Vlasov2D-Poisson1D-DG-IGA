%--------------------------------------------------------------------
% PURPOSE:
% 
% This routine computes the coefficients of the Q2 Lagrange shape 
% functions on the reference element [-1,1]x[-1,1]
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [coef]= matrix_coeff_Q2(nln)

% nln = 9

coef=zeros(nln,nln);

%   x^2    y^2   x^2y^2  x^2y   y^2x      x      y      xy      c 
M= [ 1      1      1      1      -1      -1      1      -1      1 % (-1,1)
     1      0      0      0       0      -1      0      0       1 % (-1,0)
     1      1      1      -1     -1      -1     -1      1       1 % (-1,-1)
     0      1      0      0       0       0     -1      0       1 % (0,-1)
     0      0      0      0       0       0      0      0       1 % (0,0)
     0      1      0      0       0       0      1      0       1 % (0,1)
     1      1      1      1       1       1      1      1       1 % (1,1)
     1      0      0      0       0       1      0      0       1 % (1,0)
     1      1      1      -1      1       1     -1      -1      1 % (1,-1)    
];

for i=1:nln
    f=zeros(nln,1);
    f(i)=1;
    coef(:,i)=M\f;
end

end