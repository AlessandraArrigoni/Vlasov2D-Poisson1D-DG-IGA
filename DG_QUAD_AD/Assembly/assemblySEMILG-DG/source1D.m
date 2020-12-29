function [f] = source1D(femregion, source, basis, V)
% PURPOSE : assembly of the source term of the semilagrangian DG method for
%           the test with manufactured solution for the Vlasov-Poisson problem.
%           It builds the global vector containing the source terms for the 
%           ndof_1D problems we need to solve at the same time for the different
%           velocity values.
%           For each discrete velocity value we build the source term for the 1D
%           system and then place it into the global vector.
%
% INPUT : femregion -> 1D femregion (x direction since the source term is
%           associated to the transport in x only by assumption)
%         source -> lambda function of x and y (velocity values) already
%           evaluated at the correct current time values.
%         basis -> struct containing the 1D basis functions already evaluated on
%           the chosen quadrature nodes (matrix nqn x nln = number local basis)
%         V -> struct containing the unique velocity values (from the largest positive one 
%           to the smallest negative one) and the associated 1D dof (scalar or
%           vector with 2 components)
%
% ASSUMPTION : the output vector is built according to the following order of
%              the 2D dofs: starting from top left of the domain, and moving
%              downwards row by row.

data_names = fieldnames (basis);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= basis.(data_names{iopt});']);
end

dx = w_1D*femregion.h; % The weights are already divided by 2
% Assumption: the FE space 1D is the same in both directions, so the number of
% 2D dofs is the square of the number of 1D dofs.
f = zeros(femregion.ndof^2,1);  
last_index = 0;

 
for j = 1:length(V.values)
    
    trasp = V.values(j);
    cur_dof1D = V.unique_coord_dof{j}'; % It can be either a scalar or a vector with 2 entries.
    f_sys = zeros(femregion.ndof,1); % source term of the current system for 1D problem
    
    for ie = 1:femregion.ne % loop over 1D elements
        % global indices of the basis functions on the current element
        index1D = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
        coords_elem = femregion.coords_element((ie-1)*2*ones(2,1) + [1:2]');

        pphys = coords_elem(1) + (femregion.h/2)*(nodes_1D + 1); % Coordinates of the local quadrature points 
        f_loc = zeros(femregion.nln,1);
        F = source(pphys, trasp); % vector with nqn elements
        
        for k = 1:femregion.nqn
            f_loc = f_loc + F(k)*dx(k).*dphiq(k,:)'; % column vector
        end
        
        f_sys(index1D) = f_sys(index1D) + f_loc;
    end
    
    % Build global source term
    temp = repmat(f_sys,[numel(cur_dof1D),1]); 
    f(last_index + 1 : last_index + numel(temp)) = temp;
    
    last_index = last_index + numel(temp);
end


end
