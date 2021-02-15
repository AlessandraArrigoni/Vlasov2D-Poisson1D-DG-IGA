% SOLVE_ADJOINT_LAPLACE: solve the linear system for the adjoint Laplacian problem, using hierarchical spaces.
% Assumption: same space and mesh used for the solution of the primal pde. 
% --> to be changed???
%
% The function solves the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^n)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
% p = solve_adjoint_laplace (hmsh, hspace, lap_mat, u, adjproblem_data)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the space of hierarchical splines (see hierarchical_space)
%   lap_mat: stiffness matrix assembled for the primal problem (we can use
%            it if the spaces are the same and the diffusion coefficient is
%            constant on the domain, as we assume here)
%   u: vector with the coefficients of the primal problem solution
%   adjproblem_data: a structure with data of the problem. It contains the fields:
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - f:            function handle of the source term
%    - f_type:       string denoting the source when it involves u
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - g_type:       string denoting the Neumann BC when it involves u
%    - h:            function for Dirichlet boundary condition
%
% OUTPUT:
%
%   p: computed degrees of freedom

function p = solve_adjoint_laplace (hmsh, hspace, lap_mat, u, adjproblem_data)

% Assemble source term 
if isfield(adjproblem_data, 'f')
    rhs = op_f_v_hier_adjoint(hspace, hmsh, adjproblem_data.f);
elseif (isfield(adjproblem_data, 'f_type'))
    rhs = op_f_v_hier_adjoint(hspace, hmsh, u, adjproblem_data.f_type);
else
    disp('ERROR: no source term specified for adjoint problem');
end

% Apply Neumann boundary conditions
for iside = adjproblem_data.nmnn_sides
   % Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
   if isfield(adjproblem_data, 'g')
       gside = @(varargin) adjproblem_data.g(varargin{:},iside);
       dofs = hspace.boundary(iside).dofs;
       rhs(dofs) = rhs(dofs) + op_f_v_hier_adjoint (hspace.boundary(iside), hmsh.boundary(iside), gside); 
   elseif isfield(adjproblem_data, 'g_type')
       dofs = hspace.boundary(iside).dofs;
       u_bdn = u(dofs);
       rhs(dofs) = rhs(dofs) + op_f_v_hier_adjoint (hspace.boundary(iside), hmsh.boundary(iside), u_bdn, adjproblem_data.g_type); 
   else
       disp('ERROR: no Neumann term specified for adjoint problem');
   end
end

% Apply Dirichlet boundary conditions
p = zeros (hspace.ndof, 1);
[p_dirichlet, dirichlet_dofs] = sp_drchlt_l2_proj (hspace, hmsh, adjproblem_data.h, adjproblem_data.drchlt_sides);
p(dirichlet_dofs) = p_dirichlet;

int_dofs = setdiff (1:hspace.ndof, dirichlet_dofs);
rhs(int_dofs) = rhs(int_dofs) - lap_mat(int_dofs, dirichlet_dofs)*p(dirichlet_dofs);

% Solve the linear system
p(int_dofs) = lap_mat(int_dofs, int_dofs) \ rhs(int_dofs);


end