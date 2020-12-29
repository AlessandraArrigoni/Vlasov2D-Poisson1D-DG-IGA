function [M] = global_matrixE(femregionX, femregionY, shapefun, nodes_ref, weights, dt, local_structures, E)
% PURPOSE : build the global matrix containing all the matrices (called "local"
% in the following) for the linear advection equations in v direction
% associated to the discrete values of the electric field (dofs).
% 
% We construct the matrices A0 and A1 depending on the electric field value,
% while the local structure is the same once we know the SIGN of E,
% so they are precomputed and given as input in local_structures.
%
% ASSUMPTION : we impose "COMPACT SUPPORT" boundary conditions, feature that 
%              is hardcoded in the construction of the local matrices.
%
% INPUT : E -> struct with fields "value" (array), "unique_coord_dof" (cell
%              array that associates the corresponding dof/dofs to each 
%              unique value of the transport coefficient)
%         shapefun -> contains the basis function 1D (not yet evaluated, as in
%              the output of basis_lagrange1D)
%         nodes_ref, weights -> 1D quadrature nodes and weights per element
%              in the reference element (0,1)


data_names = fieldnames (local_structures);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} ' = local_structures.(data_names{iopt});']);
end

% Vectors for the indices and values of the global matrix 
% (the local ones have all the same dimension)
rows = zeros(numel(rows_loc_minus)*femregionX.ndof, 1);
cols = zeros(numel(cols_loc_minus)*femregionX.ndof, 1);
values = zeros(numel(cols_loc_minus)*femregionX.ndof, 1);

% We must keep track of the current last nonzero entry in the global vectors
last_index = 0; 

for k = 1:length(E.values)
    
    trasp = E.values(k);
    
    % Compute values for the the local matrices (A0 and A1)
    values_loc = transport_matrix1D(femregionY, shapefun, nodes_ref, weights, trasp, dt);
    
    % Global matrix structure
    cur_dof1D = E.unique_coord_dof{k}'; % It can be either a scalar or a vector with 2 entries
    cur_values = repmat(values_loc, [numel(cur_dof1D), 1]);
    indices = reshape(repmat(cur_dof1D-1, [numel(rows_loc_minus),1]), [numel(rows_loc_minus)*numel(cur_dof1D),1]);
    
    if trasp < 0 || abs(trasp) < 1e-9
       cur_rows = indices*femregionX.ndof + repmat(rows_loc_minus, [numel(cur_dof1D),1]);
       cur_cols = indices*femregionX.ndof + repmat(cols_loc_minus, [numel(cur_dof1D),1]);
       
    elseif trasp > 0
       cur_rows = indices*femregionX.ndof + repmat(rows_loc_plus, [numel(cur_dof1D),1]);
       cur_cols = indices*femregionX.ndof + repmat(cols_loc_plus, [numel(cur_dof1D),1]);
    end
    
    rows(last_index + 1 : last_index  + numel(cur_values)) = cur_rows;
    cols(last_index + 1 : last_index  + numel(cur_values)) = cur_cols;
    values(last_index + 1 : last_index  + numel(cur_values)) = cur_values;
    
    last_index = last_index + numel(cur_values);
    
end

M = sparse(rows, cols, values);

end

