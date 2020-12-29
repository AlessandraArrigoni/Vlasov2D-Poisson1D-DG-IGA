%--------------------------------------------------------------------
% PURPOSE:
%
% This routine computes the L1 and L2 error in different norms for the linear
% transport problem in 1D, i.e.,
%
% ||u-u_h||_{1,\Omega}
% ||u-u_h||_{2,\Omega}
% 
% NB: the errors are computed using TWICE the number of quadrature nodes
% employed to compute the elements of the matrix.
%
%--------------------------------------------------------------------


function [errors]= compute_errors1D(data, femregion, uh, time)

nln = femregion.nln;
ne = femregion.ne;

% initialization
E_L2 = 0;
E_L1 = 0;
L2_norm_ex = 0;
L1_norm_ex = 0;

% scalar shape functions
[shape_basis]= basis_lagrange1D(data.fem1D);

% 1D quadrature nodes and weights: twice as many nodes.
nqn = 2*femregion.nqn;
[nodes, w] = gauleg(-1,1,nqn); % Arrays with the quadrature nodes and weights on the interval (-1,1)

% Evaluation of shape functions (they are defined on (0,1) so we need to rescale
% the gauss points
nodes_ref = (nodes+1)*0.5;
[dphiq,~] = evalshape1D(shape_basis, nodes_ref, nln); % nqn x nln

for ie = 1:ne % loop over elements

    index = (ie-1)*nln*ones(nln,1) + [1:nln]'; % indices of the 1D basis function on the current element
    coords_elem = femregion.coords_element((ie-1)*2*ones(2,1) + [1:2]'); % boundaries of the current element 

    local_uh = uh(index);
    pphys = coords_elem(1) + (femregion.h/2)*(nodes + 1); % Coordinates of the local quadrature points
    
    dx = femregion.h * w/2;
    local_exact = data.initial_f( pphys - data.trasp*time, 0); % The second parameter is y coordinate
     
    for k = 1:nqn % loop over quadrature nodes
       
        local_aprox = dphiq(k,:)*local_uh;
        
        E_L2 = E_L2 + ((local_aprox-local_exact(k)).^2).*dx(k);
        L2_norm_ex = L2_norm_ex + local_exact(k)^2.*dx(k);
        
        E_L1 = E_L1 + abs(local_aprox-local_exact(k)).*dx(k);
        L1_norm_ex = L1_norm_ex + abs(local_exact(k)).*dx(k);
                
    end

end

E_L2 = sqrt(E_L2);
L2_norm_ex = sqrt(L2_norm_ex);


errors = struct('E_L2',E_L2,...
                'L2_norm_exact', L2_norm_ex,...
                'E_L1',E_L1,...
                'L1_norm_exact', L1_norm_ex);
            
end



