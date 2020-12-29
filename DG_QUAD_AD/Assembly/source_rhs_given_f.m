function [f] = source_rhs_given_f(femregion, source, basis)
% Assembly of the rhs (source term) for a PDE with the Discontinuous Galerkin 
% method on a structured grid of triangles or quarilaterals. 
%
% INPUT: source -> lambda function f(x,v) with the expression of f in the term
%                   \int_{\Omega} f * v dx dv
%
% OUTPUT: f -> vector of length femregion.ndof with entries 
%                   \int_{\Omega} f * \phi_j dx dv
%                   where \phi_j is a DG basis function


data_names = fieldnames (basis);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= basis.(data_names{iopt});']);
end

% Get Jacobian assuming all elements have the same shape
BJ = femregion.BJ;
Jdet = abs(det(BJ));

f = sparse(femregion.ndof,1);            

for ie = 1:femregion.ne % loop over elements
    
    % Global indices of the basis functions on the current element
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';

    index_element=femregion.nedges*(ie-1).*ones(femregion.nedges,1) + [1:1:femregion.nedges]';
    
    coords_elem=femregion.coords_element(index_element, :);

    [pphys_2D] = map_local_physical_points(coords_elem, nodes_2D, BJ); 
    
    dx = w_2D*Jdet;
    
    f_loc = zeros(femregion.nln,1);
   
    F = source(pphys_2D(:,1), pphys_2D(:,2));
    
    for k = 1:length(w_2D) % loop over 2D quadrature nodes
    
        f_loc = f_loc + F(k)*dphiq(k,:)'*dx(k);
        
    end

    f(index) = f(index) + f_loc;
end


end
