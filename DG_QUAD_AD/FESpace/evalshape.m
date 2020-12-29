%--------------------------------------------------------------------
% PURPOSE:
%
% This routine evaluates the scalar shape functions and their derivatives
% at the quadrature nodes on the reference element, and on the boundary faces
% of the reference element, i.e.,
%
% dphiq  -> scalar shape functions (for volume integrals)
% B_edge -> scalar shape functions on the faces (for face integrals) CELL ARRAY
%            since its elements have different size as the number of quadrature nodes in x
%            and y can be different
% Grad   -> gradient of the scalar shape functions (for volume integrals)
% G_edge -> gradient of the scalar shape functions on the faces (for face
%           integrals). Computed only if required as output.
%
% NB : WORKS ONLY FOR CARTESIAN MESHES! not for triangle grids.
% INPUT : nodes_1D = cell array with 2 cells with the nodes on the reference
%               interval (-1,1)
%         nodes_2D = matrix nx*ny x 2 with the nodes on the reference element
%--------------------------------------------------------------------

function [dphiq, Grad, B_edge, G_edge] = evalshape(femregion,shape_basis,nodes_2D,nodes_1D,nln)

% fprintf('Evaluate the shape functions on the reference element, FEM: %s\n', femregion.fem);


nqn_x = length(nodes_1D{1});
nqn_y = length(nodes_1D{2});

c=shape_basis(1).coeff;  % needed by eval function! 

for s=1:nln % loop over the number of scalar shape functions

    csi=nodes_2D(:,1);
    eta=nodes_2D(:,2);
    dphiq(1,:,s)=eval(shape_basis(s).fbasis);
    Grad(:,1,s)=eval(shape_basis(s).Gbasis_1);
    Grad(:,2,s)=eval(shape_basis(s).Gbasis_2);

    
    % Left side
    csi = -ones(1,nqn_y);
    eta = swap_vector(nodes_1D{2});
    B_edge{1}(s,:) = eval(shape_basis(s).fbasis);

    %============================================================%
    % Bottom side
    csi = nodes_1D{1};
    eta = -ones(1,nqn_x);
    B_edge{2}(s,:) = eval(shape_basis(s).fbasis);

    %============================================================%
    % Right side
    csi = ones(1,nqn_y);
    eta = nodes_1D{2};
    B_edge{3}(s,:) = eval(shape_basis(s).fbasis);

    %============================================================%
    % Top side
    csi = swap_vector(nodes_1D{1});
    eta = ones(1,nqn_x);
    B_edge{4}(s,:) = eval(shape_basis(s).fbasis);
   
end

% Compute also G_edge only if needed.
if nargout == 4
    for s=1:nln % loop over the number of scalar shape functions
        if (femregion.nedges==3) % triangular mesh--> DOES NOT WORK
       
            % Edge 1
            csi=nodes_1D;
            eta=-zeros(1,nqn_1D);
            G_edge(:,1,s,1) = eval(shape_basis(s).Gbasis_1)';
            G_edge(:,2,s,1) = eval(shape_basis(s).Gbasis_2)';
        
            % Edge 2
            eta=nodes_1D;
            csi=ones(1,nqn_1D)-eta;
            G_edge(:,1,s,2) = eval(shape_basis(s).Gbasis_1)';
            G_edge(:,2,s,2) = eval(shape_basis(s).Gbasis_2)';

            % Edge 3
            csi=zeros(1,nqn_1D);
            eta=ones(1,nqn_1D)-nodes_1D;
            G_edge(:,1,s,3) = eval(shape_basis(s).Gbasis_1)';
            G_edge(:,2,s,3) = eval(shape_basis(s).Gbasis_2)';

        elseif (femregion.nedges==4) % quad mesh
            
            % Left side
            csi = -ones(1,nqn_y);
            eta = swap_vector(nodes_1D{2});
            G_edge(:,1,s,1) = eval(shape_basis(s).Gbasis_1)';
            G_edge(:,2,s,1) = eval(shape_basis(s).Gbasis_2)';
            
            % Bottom side
            csi = nodes_1D{1};
            eta = -ones(1,nqn_x);
            G_edge(:,1,s,2) = eval(shape_basis(s).Gbasis_1)';
            G_edge(:,2,s,2) = eval(shape_basis(s).Gbasis_2)';

            % Right side
            csi = ones(1,nqn_y);
            eta = nodes_1D{2};
            G_edge(:,1,s,3) = eval(shape_basis(s).Gbasis_1)';
            G_edge(:,2,s,3) = eval(shape_basis(s).Gbasis_2)';

            % Top side
            csi = swap_vector(nodes_1D{1});
            eta = ones(1,nqn_x);
            G_edge(:,1,s,4) = eval(shape_basis(s).Gbasis_1)';
            G_edge(:,2,s,4) = eval(shape_basis(s).Gbasis_2)';
    
        end
    end
end


% fprintf('--- END --- Evaluate the shape functions on the reference element, FEM: %s\n', femregion.fem);
end