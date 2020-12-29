function [f] = rhs_givenf1D(femregion, source, basis)
%--------------------------------------------------------------------
% PURPOSE:
%
% Assembly of the right hand side for the discontinuous galerkin method 1D. 
%
% INPUT: basis -> basis functions 1D evaluated at quadrature points on the 
%                 reference element (0,1)
%        source -> lambda function for the source term
%
%--------------------------------------------------------------------------


%fprintf('--------Begin computing rhs 1D for direction %s\n', femregion.direction);

data_names = fieldnames (basis);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} ' = basis.(data_names{iopt});']);
end

nln = femregion.nln;
dx = w_1D*femregion.h; % The weights are already divided by 2
f = zeros(femregion.ndof,1);

for ie = 1:femregion.ne % loop over elements

    index = (ie-1)*nln*ones(nln,1) + [1:nln]';
    coords_elem = femregion.coords_element((ie-1)*2*ones(2,1) + [1:2]');

    pphys = coords_elem(1) + (femregion.h/2)*(nodes_1D + 1); % Coordinates of the local quadrature points;
    F = source( pphys ); 
     
    f_loc = zeros(nln,1);
    
    for k = 1:length(nodes_1D) % loop over quadrature nodes
        f_loc = f_loc + F(k)*dphiq(k,:)'.*dx(k);                
    end
    
    f(index) = f(index) + f_loc;

end

end