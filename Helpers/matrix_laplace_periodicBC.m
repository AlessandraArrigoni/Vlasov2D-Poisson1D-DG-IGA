function [geometry, msh, space, lap_mat] = matrix_laplace_periodicBC (problem_data, method_data)
% Returns the matrix lap_mat and the geometry/space structures of GeoPDEs 
% used to solve the 1D poisson problem:
% - div ( epsilon(x) grad (u)) = f   
% in Omega = F((0,1)^n) 
% with periodic boundary conditions (user defined regularity) and
% imposing an average on the solution via a Lagrange multiplier (if required:
% tested with average = 0 only)


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

% Perform degree elevation and knot insertion for the unclamping (needed for
% periodic spaces)
domain_elev = nrbdegelev(geometry.nurbs, degree - 1); % number of times we need to increase the order (assuming we start from linear)
[knots, zeta] = kntrefine (domain_elev.knots, nsub-1, degree, regularity);
domain_clamped = nrbkntins(domain_elev, knots(degree+2 : end-degree-1));
domain_unclamped = nrbunclamp(domain_clamped, continuityBC);
  
% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

% Construct space structure: the space has the same number of basis functions as
% if it were not periodic. Only the matrix dimensions are affected by the
% reduction of the basis.
space    = sp_bspline (domain_unclamped.knots, degree, msh);

% Compute matrices' values
[rows, cols, values] = op_gradu_gradv_tp (space, space, msh, c_diff);

% Apply periodic boundary conditions
% OLD VERSION
% for j = 0 : continuityBC
%     jj = continuityBC - j;
%     rows(rows == space.ndof - jj) = 1 + j;
%     cols(cols == space.ndof - jj) = 1 + j;
%    
% end

% NEW VERSION (FASTER)
indices_new = [1 : space.ndof - continuityBC - 1, 1 : continuityBC + 1];
rows = indices_new(rows);
cols = indices_new(cols);

% Assembly the matrices: the last two parameters of function sparse is the
% number of rows and columns, i.e. the number of independent basis functions,
% equals to ndof - k - 1 if we impose continuity k at the boundaries. 

if set_avg
    % The solution is constrained to have 0 average via a Lagrange multiplier
    % that adds a row and a column to the system
    B = op_f_v_tp(space, msh, @(x) 1 + 0.*x); % column vector
        
    % "Apply periodic BC" also to B since the basis functions we consider are the periodic ones
    % OLD VERSION
    %     for j = 0 : continuityBC
    %         jj = continuityBC - j;
    %         B(1 + j)= B(1 + j) + B(space.ndof - jj);    
    %     end
    
    % NEW VERSION (FASTER)
    B(1 : continuityBC + 1) = B(1 : continuityBC + 1) + B(end - continuityBC : end);
    
    B(space.ndof - continuityBC : space.ndof) = []; % remove additional entries 
    
    lap_mat = sparse(rows, cols, values, space.ndof - continuityBC, space.ndof - continuityBC);
    lap_mat(1:end-1,end) = B;
    lap_mat(end, 1:end-1) = B';
    
else
    % In this case the solution is not constrained
    lap_mat = sparse(rows, cols, values, space.ndof - continuityBC - 1, space.ndof - continuityBC - 1);
    
end


end