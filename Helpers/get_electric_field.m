% PURPOSE: Compute the electric field starting from the solution to Poisson
% problem \phi and taking its derivative.
% ASSUMPTION:  E = \grad \phi -> they have the same sign.
%
% INPUT: space, msh -> structures from geoPDEs. Run precompute function before!
%        u_poisson -> vector with the coeffiecient of the solution to Poisson pb
% OUTPUT: ef -> values of the electric field, matrix n_quad_nodes_per_el X n_el
%         xpoints -> coordinates of the 1D quadrature nodes.

function [ef, xpoints] = get_electric_field(space, msh, u_poisson)

ef = zeros(msh.nqn_dir(1), msh.nel_dir(1));

uc = zeros (size (space.connectivity)); % nqn x nel
uc(space.connectivity~=0) = ...
        u_poisson(space.connectivity(space.connectivity~=0));
    
for iel = 1:msh.nel_dir(1)
    basis = permute(space.shape_function_gradients(:,:,:,iel), [2,3,1]); % Shape functions evaluated on each element
    ef(:,iel) = basis * uc(:,iel);
end
xpoints = permute(msh.geo_map, [2,3,1]);


end
