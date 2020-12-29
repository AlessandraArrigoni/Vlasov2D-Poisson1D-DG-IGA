%--------------------------------------------------------------------
% PURPOSE:
%
% This routine computes the errors on the DG solution in L2 norm, i.e.,
%
% ||u-u_h||_{0,\Omega}
% 
% NB: the errors are computed using TWICE as many quadrature nodes
% employed to compute the elements of the matrix 
% 
% INPUT: solution -> discrete solution of Vlasov equation: vector with the
%           coefficients with respect to the DG basis (femregion.ndof)
%
% OUTPUT: errors -> struct containing the global L2 norm of the error and the L2
%           norm of the exact solution.
%
%--------------------------------------------------------------------


function [errors]= compute_errors(Data, femregion, solution)

nln=femregion.nln;
ne=femregion.ne;

% Initialization
E_L2 = 0;
L2_norm_ex = 0;

% Get Jacobian: assuming all the elements have the same shape
BJ = femregion.BJ;
Jdet = abs(det(BJ));

% Scalar shape functions
[shape_basis]= basis_lagrange(Data.fem);

% 1D and 2D quadrature nodes and weights 
nqn = 2*Data.nqn;
[nodes_1D, ~, nodes_2D, w_2D] = quadrature(nqn);
% nodes_1D and w_1D are cell arrays

% Evaluation of shape functions
[dphiq, ~, ~]= evalshape(femregion, shape_basis, nodes_2D, nodes_1D, femregion.nln);
dphiq = reshape(dphiq,[length(nodes_2D),femregion.nln]); % nqn2D x nln

dx = Jdet*w_2D;

for ie=1:ne % loop over elements

    % global indices of the basis functions on the current element
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    index_element=femregion.nedges*(ie-1).*ones(femregion.nedges,1) + [1:1:femregion.nedges]';

    coords_elem = femregion.coords_element(index_element, :);

    % Local solution vector (coefficients)
    local_uh = solution(index);
    
    % Jacobian and physical coordinates at the quadrature points
    [pphys_2D] = map_local_physical_points(coords_elem, nodes_2D, BJ); 
    
    if isfield(Data,'test') && strcmp(Data.test, 'solid_body')
        rvett = sqrt((pphys_2D(:,1)-Data.x0).^2+(pphys_2D(:,2)-Data.y0).^2);
        local_exact = Data.initial_f(rvett, 0); 
        % It only works if we compute the error at a point where the solution is
        % equal to the initial one.
    else
        t = Data.time;
        local_exact = Data.exact_sol(pphys_2D(:,1), pphys_2D(:,2),t);
    end
    
    for k = 1:length(w_2D) % loop over quadrature nodes

        local_aprox = dphiq(k,:)*local_uh;
        
        E_L2 = E_L2 + ((local_aprox - local_exact(k)).^2).*dx(k);
        L2_norm_ex = L2_norm_ex + local_exact(k)^2.*dx(k);
                
    end

end

E_L2 = sqrt(E_L2);
L2_norm_ex = sqrt(L2_norm_ex);

errors = struct('E_L2',E_L2,...
    'L2_norm_exact', L2_norm_ex);

end



