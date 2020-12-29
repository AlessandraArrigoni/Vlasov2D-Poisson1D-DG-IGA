function [rows, cols, values, valuesinv] = mass_matrix1D_structure(femregion, basis)
%--------------------------------------------------------------------
% PURPOSE:
%
% Assembly of the MASS matrix M and its inverse for the discontinuous Galerkin 
% method 1D (block diagonal). The values and the corresponding indices of rows
% and columns are returned in vectors to be used by sparse().
%
% Since the entries of the local M and Minv are the same for all elements, we
% compute them only once and then build the global matrix by positioning them in
% the right place.
%
% INPUT: basis -> contains the basis functions 1D evaluated at quadrature
%                 points on the reference element [0,1]
%--------------------------------------------------------------------------


% Read fields of the inputs into local variables
data_names = fieldnames (basis);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} ' = basis.(data_names{iopt});']);
end

% Useful constants to build the vectors with the rows and columns' indices
% needed to assemble the sparse matrix at the end of the function.
nln = femregion.nln;
nmax_el = nln^2; % nln columns * nln rows
nmax = femregion.ne*nmax_el ;
rows = zeros(nmax, 1);
cols = zeros(nmax, 1);
values = zeros(nmax, 1);
valuesinv = zeros(nmax, 1);

dx = w_1D*femregion.h; % The weights are already divided by 2

M_loc = zeros(femregion.nln, femregion.nln);
for k = 1:length(dx)
    M_loc = M_loc + dphiq(k,:)'*dphiq(k,:).*dx(k);
end
values_loc = reshape(M_loc, [nmax_el, 1]);

% Compute local inverse
Minv_loc = M_loc \ eye(nln);
valuesinv_loc = reshape(Minv_loc, [nmax_el,1]);

for ie=1:femregion.ne % loop over elements
    
    % Global indices of the basis functions on the current element
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
 
    % Assemble vectors for global matrices
    values((ie-1)*nmax_el + 1 : ie*nmax_el) = values_loc;
    
    rows_loc = repmat(index, [nln,1]);
    rows((ie-1)*nmax_el + 1 : ie*nmax_el) = rows_loc;
    
    cols_loc = reshape(ones(nln,1)*index', [nmax_el, 1]);
    cols((ie-1)*nmax_el + 1 : ie*nmax_el) = cols_loc;
    
    valuesinv((ie-1)*nmax_el + 1 : ie*nmax_el) = valuesinv_loc;
       
end


end