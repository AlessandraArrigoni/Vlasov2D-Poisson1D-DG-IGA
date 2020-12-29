%%--------------------------------------------------------------------
% PURPOSE:
%
% This routine computes the [nx, ny] Gauss-Ledendre nodes and weights on the
% reference interval (-1,1)  to be used for face integrals, and
% the nx*ny Gauss-Legendre nodes and weights on the reference square
% (-1,1)x(-1,1) to be used for volume integrals.
%
% INPUT: n = [nx, ny] vector specifying the number of quadrature nodes in the
% first and second direction of the domain.
%
% OUTPUT: node_1D = cell array with node_1D{1} = nodes in x, node_1D{2} = nodes in y
%         w_1D = cell array with the weights already divided by 2 (rescaling)
%         node_2D = matrix nx*ny x 2 (for each x, we put all the ys then change x)
%         w_2D = vector nx*ny
%--------------------------------------------------------------------

function [node_1D, w_1D, node_2D, w_2D] = quadrature(n)

        [nodesx, wx]= gauleg(-1,1,n(1));
        [nodesy, wy]= gauleg(-1,1,n(2));
        
        node_1D{1} = nodesx;
        node_1D{2} = nodesy;
        
        w_1D{1} = wx./2;
        w_1D{2} = wy./2;
                
        node_2D = [];
        w_2D = [];
        for i=1:n(1)
            for j=1:n(2)
                node_2D = [node_2D; nodesx(i), nodesy(j)];
                w_2D = [w_2D, wx(i).*wy(j)];
            end
        end
end


