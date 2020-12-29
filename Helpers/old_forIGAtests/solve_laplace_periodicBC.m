% SOLVE_LAPLACE: Solve a Laplace problem with a B-spline discretization (non-isoparametric approach). 
%
% The function solves the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^n)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_laplace_periodicBC (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - f:            source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete space (see sp_scalar)
%  u:        the computed degrees of freedom
%
% See also EX_LAPLACE_SQUARE, EX_LAPLACE_THICK_RING for examples.
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, 2015, 2017 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [geometry, msh, space, u] = ...
              solve_laplace_periodicBC (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure
geometry  = geo_load (geo_name);

% Perform degree elevation and knot insertion for the unclamping
domain_elev = nrbdegelev(geometry.nurbs, degree - 1); % how many times we need to increase the order (assuming we start from linear)
[knots, zeta] = kntrefine (domain_elev.knots, nsub-1, degree, regularity);
domain_clamped = nrbkntins(domain_elev, knots(degree+2 : end-degree-1));
domain_unclamped = nrbunclamp(domain_clamped, continuityBC);
  
% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

% Construct space structure
space    = sp_bspline (domain_unclamped.knots, degree, msh);
space_cla = sp_bspline(domain_clamped.knots, degree, msh);

% Plot basis functions
if(plot_basis_fun)
    pts = linspace(0,1,200);
    figure(1)
    for i  = 1:space.ndof
        cfs = zeros(1, space.ndof);
        cfs(i) = 1;    
        [basis, F] = sp_eval(cfs, space, geometry, {pts});
        figure(1)
        plot(F, basis, '-','LineWidth',1)
        hold on
        grid on
       
    end
    
    figure(2)
    for i  = 1:space.ndof
        cfs = zeros(1, space.ndof);
        cfs(i) = 1;    
        [basis, F] = sp_eval(cfs, space_cla, geometry, {pts});
        figure(2)
        plot(F, basis, '-','LineWidth',1)   
        hold on
        grid on
    end
end

% Compute matrices' values
[rows, cols, values] = op_gradu_gradv_tp (space, space, msh, c_diff);
if isa(f, 'function_handle')
    rhs       = op_f_v_tp (space, msh, f); 
else
    rhs = op_f_v_tp_data(space, msh, f); % se la f è la matrice calcolata come integrale di linea
end

% Apply periodic boundary conditions
% Se accedo a un vettore "corto" con un vettore di indici più lungo, ottengo 
% un vettore "lungo" con valori ripetuti in base a quanto detto dagli indici
indices_new = [1 : space.ndof - continuityBC - 1, 1 : continuityBC + 1];
rows = indices_new(rows);
cols = indices_new(cols);

rhs(1 : continuityBC + 1) = rhs(1 : continuityBC + 1) + rhs(end - continuityBC : end);
% VERSIONE VECCHIA 
% for j = 0 : continuityBC
%     jj = continuityBC - j;
%     rows(rows == space.ndof - jj) = 1 + j;
%     cols(cols == space.ndof - jj) = 1 + j;
%     
%     rhs(1 + j)= rhs(1 + j) + rhs(space.ndof - jj);    
% end


% Assembly the matrices: the last two parameters of function sparse is the
% number of rows and columns, i.e. the number of independent basis functions,
% equals to ndof - k - 1 if we impose continuity k at the boundaries. 
if set_avg
    % The solution is constrained to have 0 average via a Lagrange multiplier
    % that adds a row and a column to the system
    B = op_f_v_tp(space, msh, @(x) 1 + 0.*x); % column vector
     
    % "Apply periodic BC" also to B since the basis functions we consider are the periodic ones
    B(1 : continuityBC + 1) = B(1 : continuityBC + 1) + B(end - continuityBC : end);
    % Versione vecchia
%     for j = 0 : continuityBC
%         jj = continuityBC - j;
%         B(1 + j)= B(1 + j) + B(space.ndof - jj);    
%     end
    
    B(space.ndof - continuityBC : space.ndof) = [];
    stiff_mat = sparse(rows, cols, values, space.ndof - continuityBC, space.ndof - continuityBC);
    stiff_mat(1:end-1,end) = B;
    stiff_mat(end, 1:end-1) = B';
    
    rhs(space.ndof - continuityBC + 1 : space.ndof) = []; % reduce the rhs 
    rhs(end) = avg; % equation associated to the zero average constraint
    
    % Solve the linear system
    u = stiff_mat \ rhs;
    u(space.ndof - continuityBC : space.ndof) = u(1:continuityBC + 1); % add the values associated to the dofs on the other side of the domain
   
    % NB: avremmo dovuto togliere l'ultimo elemento della soluzione (quello
    % associato al moltiplicatore di Lagrange) e poi aggiungere i valori secondo
    % le BC periodiche, però provo a fare tutto in un colpo, sperando di non
    % sbagliare indici.
    
else
    % In this case the solution is not constrained
    stiff_mat = sparse(rows, cols, values, space.ndof - continuityBC - 1, space.ndof - continuityBC - 1);
    rhs(space.ndof - continuityBC : space.ndof) = []; % reduce the rhs

    % Solve the linear system
    u = stiff_mat\rhs;
    u(space.ndof - continuityBC : space.ndof) = u(1:continuityBC + 1); % add the values associated to the dofs on the other side of the domain
end



end
