function [M, Minv] = mass_matrix(femregion, basis)
%--------------------------------------------------------------------
% PURPOSE:
%
% Assembly of the MASS matrix M for the discontinuous galerkin method (block diagonal).
% 
% No boundary condition is imposed to the matrix since it involves only volume
% integrals and does not interact with the neighbours.
%
% OUTPUT: M -> mass matrix
%         Minv -> inverse of mass matrix (computed only if required as output)
%--------------------------------------------------------------------


fprintf('--------Begin computing MASS matrix for Q%d\n',femregion.degree);

data_names = fieldnames (basis);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} ' = basis.(data_names{iopt});']);
end

% Useful constants to build the vectors with the rows and columns' indices
% needed to assemble the sparse matrix at the end of the function.
nln = femregion.nln;
nmax_el = nln^2; % nln colonne * nln righe -> max number of non zero entries for each element
nmax = femregion.ne*nmax_el ;
rows = zeros(nmax, 1);
cols = zeros(nmax, 1);
values = zeros(nmax, 1);
valuesinv = zeros(nmax, 1);


% Get Jacobian: assuming all the elements of the mesh have the same shape.
BJ = femregion.BJ;
Jdet = abs(det(BJ));

for ie=1:femregion.ne % loop over elements
    
    % Global indices of the basis functions on the current element
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
      
    dx = w_2D*Jdet;
    
    M_loc = zeros(femregion.nln, femregion.nln);
    
    for k=1:length(w_2D) % loop over 2D quadrature nodes
        M_loc = M_loc + dphiq(k,:)'*dphiq(k,:).*dx(k);        
    end
    
    % Assemble global vectors 
    values_loc = reshape(M_loc, [nmax_el,1]); % vector with values
    values((ie-1)*nmax_el + 1 : ie*nmax_el) = values_loc;
    
    rows_loc = repmat(index, [nln,1]);
    rows((ie-1)*nmax_el + 1 : ie*nmax_el) = rows_loc;
    
    cols_loc = reshape(ones(nln,1)*index', [nmax_el, 1]);
    cols((ie-1)*nmax_el + 1 : ie*nmax_el) = cols_loc;
    
    if nargout == 2 
        Minv_loc = M_loc \ eye(length(index));
        valuesinv_loc = reshape(Minv_loc, [nmax_el,1]);
        valuesinv((ie-1)*nmax_el + 1 : ie*nmax_el) = valuesinv_loc;
    end
    
    
end

% Assemble global matrices from the vectors
M = sparse(rows, cols, values);
if nargout == 2
    Minv = sparse(rows, cols, valuesinv);
end

fprintf('--------End computing MASS matrix for Q%d\n', femregion.degree);

end