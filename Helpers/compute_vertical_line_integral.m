function [values, nodes_glob] = compute_vertical_line_integral(basis, u_h, femregion, region)
% Compute line integral of u_h along the "vertical": \int_{\Omega_v} f(x,v,t) dv
% where f is the numerical solution to Vlasov equation (in our case) 
% INPUT : 
%         basis -> struct with DG basis functions evaluated on the 2D quadrature nodes
%         u_h -> coefficient of the solution f with respect to the DG basis
%           (vector with femregion.ndof elements) 
%         femregion 
%         region -> contains information on the 2D elements sharing the same
%           x-interval (see region.column_elements)
%
% OUTPUT : 
%         values -> matrix with the values of the line integral for each discrete x
%           (that is, for each 1D quadrature point in the x direction)
%            n_rows = n_nodes1D per element, n_columns = n_elements in x dir.
%         nodes_glob -> vector with the 1D quadrature nodes in order


data_names = fieldnames (basis);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= basis.(data_names{iopt});']);
end

% u_h must be a column vector
if size(u_h,2) > size(u_h,1)
    u_h = u_h';
end

nqn_x = length(w_1D{1});
nqn_y = length(w_1D{2});
values = [];
nodes_glob = [];

% Loop on 1D elements
for elx = 1:numel(region.x_points1D)-1
    
    % Create local vector to store the integrals
    values_loc = zeros(nqn_x,1);
    
    % Create local vector with the 1D quadrature nodes (used for Poisson pb.)
    hh = region.x_points1D(elx+1) - region.x_points1D(elx);
    node1D_loc = region.x_points1D(elx) +0.5*hh*(nodes_1D{1}+1); % translate nodes from [-1,1] to current interval
    
    % Loop on the associated 2D elements
    for ie = region.column_elements(elx,:)
        
        % Find 2D basis functions living on the current element
        idx_shape = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
        u_h_loc = u_h(idx_shape);
        
        % Find vertical length of the current element (assumption: rectangle)
        index_element = femregion.nedges*(ie-1).*ones(femregion.nedges,1) + [1:1:femregion.nedges]';
        coords_elem = femregion.coords_element(index_element, :);
        hy = max(coords_elem(:,2)) - min(coords_elem(:,2)); % ymax - ymin
        
        % Compute function f on the 2D quadrature nodes (column vector)
        sol_loc = dphiq*u_h_loc;
        
        % Compute line integral with 1D weights (already divided by 2, so we
        % just need the element length to scale them correctly)
        for k = 1:nqn_x
            values_loc(k) = values_loc(k) + (w_1D{2}*sol_loc((k-1)*nqn_y + 1 : k*nqn_y))*hy;
        end
        
    end
    
    % Add to the output matrix
    nodes_glob = [nodes_glob, node1D_loc];
    values = [values, values_loc];
end

end

