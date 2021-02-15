function [dofs] = select_dofs_from_shapeder(dJ, hspace, thresh)

% SELECT_DOFS_FROM_SHAPEDER : Select the dofs to keep in the refined geometric 
% space acccording to the value of the associated shape derivative.
%
% INPUT:
%
%   dJ:    vector hspace.ndof x 2 collecting the values of the shape derivatives in the
%          direction of all the dofs in hspace, componentwise. It contains NaN in
%          the indices of the bulk dofs.
%   hspace:  hierarchical space for the geometry obtained by dyadic
%            refinement on the boundary only.
%   thresh: scalar 
%
% OUTPUT:
%
%   dofs: cell array with the indices of the needed dofs per level, or []


dofs = cell(hspace.nlevels, 1);
last_idx = 0;

% Global indices of the larger shape derivatives
% We only consider the second component for the moment (sin example)
large_dJ = find(dJ(:,2) > thresh); 
% indices_nan = find(isnan(dJ(:,2))); 

for lev = 1:hspace.nlevels
    [~, ~, dofs_boundary] = intersect(large_dJ, [last_idx + 1 : last_idx + hspace.ndof_per_level(lev)]);
%     [~, ~, dofs_bulk] = intersect(indices_nan, [last_idx + 1 : last_idx + hspace.ndof_per_level(lev)]);
     
    % Indices in the tensor product spaces
    dofs_boundary_tp = hspace.active{lev}(dofs_boundary);
%     dofs_bulk_tp = hspace.active{lev}(dofs_bulk);
%     dofs{lev} = union(dofs_boundary_tp, dofs_bulk_tp);
    dofs{lev} = dofs_boundary_tp;
    last_idx = last_idx + hspace.ndof_per_level(lev);
end

end

