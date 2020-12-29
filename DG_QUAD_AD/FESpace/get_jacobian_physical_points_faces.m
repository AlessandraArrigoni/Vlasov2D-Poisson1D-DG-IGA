% For each face compute the jacobian of the mapping from the refernce
% points to the physical point
% INPUT : node_1D = cell array with 2 cells: nodes along x and nodes along y
%
% OUTPUT : pphys_1D = cell array with 4 cells (one per face): each one is a
%       matrix nx x 2 or ny x 2 with the physical coordinates of the quadrature nodes on this face. 

function [BJ_face, BJ_face_inv, pphys_1D] = get_jacobian_physical_points_faces(loc_coord, node_1D, type_mesh)

nqn_x = length(node_1D{1});
nqn_y = length(node_1D{2});
nfaces = length(loc_coord(:,1));

% Only valid for CARTESIAN meshes
x0=loc_coord(1,1);   % x-coordinates of vertices
x1=loc_coord(2,1);
x2=loc_coord(3,1);
x3=loc_coord(4,1);

y0=loc_coord(1,2);   % y-coordinates of vertices
y1=loc_coord(2,2);
y2=loc_coord(3,2);
y3=loc_coord(4,2);

% Old version without cell arrays
% nodes_face(:,:,1)=[-ones(1,nqn_y);swap_vector(node_1D{2})]';
% nodes_face(:,:,2)=[node_1D{1};-ones(1,nqn_x)]';
% nodes_face(:,:,3)=[ones(1,nqn_y);node_1D{2}]';
% nodes_face(:,:,4)=[swap_vector(node_1D{1});ones(1,nqn_x)]';

% Numbering of the faces: [1,2,3,4] = [left, bottom, right, top]
nodes_face{1} = [-ones(1,nqn_y);swap_vector(node_1D{2})]'; 
nodes_face{2} = [node_1D{1};-ones(1,nqn_x)]';
nodes_face{3} = [ones(1,nqn_y);node_1D{2}]';
nodes_face{4} = [swap_vector(node_1D{1});ones(1,nqn_x)]';

% Helper structure: associates to each face the corresponding number of nodes
n_nodes = [nqn_y, nqn_x, nqn_y, nqn_x];

for iedg = 1:nfaces
       
    BJ_face(:,:,iedg)  = (0.25) .* [-x0 - x1 + x2 + x3 ,  x0 - x1 - x2 + x3 ; -y0 - y1 + y2 + y3 , y0 - y1 - y2 + y3];
    trans = (0.25) .* [ x0 + x1 + x2 + x3 ;  y0 + y1 + y2 + y3];                       % translation vector
    BJ_face_inv(:,:,iedg) = inv(BJ_face(:,:,iedg)); % inverse of Jacobian of elemental map

   
    for k = 1:n_nodes(iedg)
        pphys_1D{iedg}(k,:) = transpose((BJ_face(:,:,iedg)*transpose(nodes_face{iedg}(k,:))+trans));
    end

end


