function [u] = solve_laplace1D_periodicBC_given_matrix(f, L, U, space, msh, method_data, problem_data)
% Solves the 1D Poisson prolem:
% - div ( epsilon(x) grad (u)) = f    
% in Omega = F((0,1)^n) 
% with periodic boundary conditions (user defined regularity) and
% imposing 0 average on the solution via a Lagrange multiplier (if required).
%
% The matrix is given as input as well as the geometry / space structures. 
% L, U are in fact the lower and upper triangular matrix obtained by LU
% factorization.
%
% The SOURCE TERM f can be either a function handle, or a matrix 
% (nrows=nqn per element, ncols=nelements 1D) with the values on the quadrature nodes.


% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Build source term
if isa(f, 'function_handle')
    rhs = op_f_v_tp (space, msh, f); 
else
    rhs = op_f_v(space, msh, f); % if f is a matrix with the values on the quadrature nodes
end

% Apply periodic boundary conditions
% OLD VERSION
% for j = 0 : continuityBC
%     jj = continuityBC - j;
%     
%     rhs(1 + j)= rhs(1 + j) + rhs(space.ndof - jj);    
% end

% NEW VERSION (FASTER)
rhs(1 : continuityBC + 1) = rhs(1 : continuityBC + 1) + rhs(end - continuityBC : end);

% Impose given average on the solution
if set_avg
    rhs(space.ndof - continuityBC + 1 : space.ndof) = []; % reduce the rhs 
    rhs(end) = avg; % equation associated to the zero average constraint
    
    % Solve the linear system: L is lower triangular such that LU = lap_mat
    y = L\rhs;
    u = U\y;
    % Add the values associated to the dofs on the other side of the domain
    % (the space contains the same number of basis functions as if it were not
    % periodic, so we associate the same coefficient u to the functions that are
    % actually the same in the periodic setting)
    u(space.ndof - continuityBC : space.ndof) = u(1:continuityBC + 1); 
   
else
    rhs(space.ndof - continuityBC : space.ndof) = []; % reduce the rhs

    % Solve the linear system
    u = stiff_mat\rhs;
    % add the values associated to the dofs on the other side of the domain
    u(space.ndof - continuityBC : space.ndof) = u(1:continuityBC + 1); 

end

end