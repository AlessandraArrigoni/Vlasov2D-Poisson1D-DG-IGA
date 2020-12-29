function [values] = transport_matrix1D(femregion1D, basis1D, nodes1D, w1D, trasp, dt)
% ---------------------------------------------------------------------
% PURPOSE:
% Compute the values of global matrix needed to make the solution advance 
% from t^n to t^(n+1), given the value of the velocity (a) of the 1D linear
% transport problem. (see. Crouseilles pg 215).
%
% STEPS:
% 1) compute parameter \alpha under the assumption that a = 0 gives \alpha = 0. 
% 2) evaluate the 1D basis functions at the correct quadrature points.
% 3) compute the local matrices A0 and A1.
% 4) create the output vector with the values of the global circulant matrix
%    in the correct order according to the sign of the transport coefficient.
%
% STRONG ASSUMPTION: in the X-direction we always impose periodic boundary
% conditions, while in the Y-direction we always impose "compact support"
% boundary conditions. This is the case for the considered tests for
% Vlasov-Poisson system, so it is hardcoded, but prevents the generalization of
% the code as it is to other problems.
% 
% INPUT : 
%   basis1D -> struct containing the basis functions 1D associated to the chosen
%              polynomial space
%   nodes1D -> array containing the quadrature nodes on the reference interval (-1,1)
%   w1D -> array containing the quadrature weights on the interval (0,1)
%   trasp -> value of the transport coefficient
%   dt -> time discretization parameter (chosen in order to satisfy the
%         assumption that the characteristics' starting point must be in the
%         neighbouring cell at most)
% ---------------------------------------------------------------------

% Compute parameter \alpha
if abs(trasp) <= 1e-9 % if trasp is small, set it to 0
    trasp = 0;
end
if trasp <= 0
    alpha = -trasp*dt/femregion1D.h;
elseif trasp > 0
    alpha = 1 - trasp*dt/femregion1D.h;
end
if (alpha < 0 || alpha > 1)
    error('Parameter alpha must be in [0,1)! Reduce dt!');
elseif (abs(alpha - 1) < 1e-11)
    error('Parameter alpha must be smaller than 1!');
end

% Evaluate shape functions: test and trial must be evaluated at different points
% depending on A0 and A1.
nodes_trialA0 = alpha + (1-alpha).*nodes1D;
[basis_trialA0, ~] = evalshape1D(basis1D, nodes_trialA0, femregion1D.nln);

nodes_testA0 = (1-alpha).*nodes1D;
[basis_testA0, ~] = evalshape1D(basis1D, nodes_testA0, femregion1D.nln);

nodes_trialA1 = alpha.*nodes1D;
[basis_trialA1, ~] = evalshape1D(basis1D, nodes_trialA1, femregion1D.nln);

nodes_testA1 = alpha.*(nodes1D - 1) + 1;
[basis_testA1, ~] = evalshape1D(basis1D, nodes_testA1, femregion1D.nln);

% Compute local values for the matrices A0 and A1
A0 = zeros(femregion1D.nln, femregion1D.nln);
for k = 1:length(w1D)
    A0 = A0 + basis_testA0(k,:)'*basis_trialA0(k,:).*w1D(k); % integrals on (0,1)
end
A0 = A0*femregion1D.h*(1-alpha); % rescaling due to the change of variable

A1 = zeros(femregion1D.nln, femregion1D.nln);
for k = 1:length(w1D)
    A1 = A1 + basis_testA1(k,:)'*basis_trialA1(k,:).*w1D(k); % integrals on (0,1)
end
A1 = A1*femregion1D.h*alpha; % rescaling due to the change of variable

% Create vector with the values for the global circulant matrix according to the
% different boundary conditions.
nln = femregion1D.nln;
nmax_el = 2*nln^2; % 2*nln columns * nln rows

if femregion1D.direction == 'x' % Periodic BC
    values_loc = reshape([A0, A1], [nmax_el,1]); % reshape takes the values columnwise
    values = repmat(values_loc, [femregion1D.ne,1]);
    
elseif femregion1D.direction == 'y' % "Compact support" BC
    values_loc = reshape([A0, A1], [nmax_el,1]); % reshape takes the values columnwise
    values = repmat(values_loc, [femregion1D.ne - 1,1]); 
    values_A0 = reshape(A0, [numel(A0), 1]);
    values_A1 = reshape(A1, [numel(A1), 1]);
    
    if trasp <= 0
        % complete values vector
        values = [values; values_A0];
    elseif trasp > 0
        values = [values_A1; values];
    end
    
end


end

