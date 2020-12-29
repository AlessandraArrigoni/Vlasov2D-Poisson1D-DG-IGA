%--------------------------------------------------------------------
% PURPOSE:
%
% This routine generates a structure containing all the information of the
% DG finite element space, i.e.,
%
% femregion.fem -> finite element space
% femregion.domain -> coordinates of the (square) domain
% femregion.type_mesh -> type of mesh
% femregion.h -> mesh size
% femregion.nln -> number of local nodes
% femregion.ndof -> total number of degrees of freedom
% femregion.ne -> total number of elements
% femregion.dof -> coordinates of the degrees of freedom
% femregion.nqn -> number of 1D quadrature nodes
% femregion.degree -> polynomial degree
% femregion.coord -> coordinates of the nodes of the mesh
% femregion.connectivity -> connectivity matrix for the mesh
% femregion.coords_element -> coordinates of the P1-DG degrees of freedom
% femregion.boundary_edges -> boundary edges
% femregion.column_elements -> matrix with list of elements that share the same
%                           x-interval
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [femregion]=create_dof(Dati,region)

fprintf('Computing dof for elements %s\n',Dati.fem);

nln=0;
ne=region.ne;
nedge=region.nedge;

dof_hat =[];

switch Dati.fem  % dof on the reference element

    case{'P1'}
        degree=1;
        dof_hat =[ 0 1 0
            0 0 1]';
        
    case{'P2'}
        degree=2;
        dof_hat =[ 0    0.5     1   0.5     0   0
                   0    0       0   0.5     1   0.5]';
               
    case{'P3'}
        degree=3;
        dof_hat =[ 0   0.25   0.75  1   0.75  0.25  0  0    0    1/3
                   0   0      0     0   0.25  0.75  1  0.75 0.25 1/3]';   

    case{'Q1'}
        degree=1;
        dof_hat=[-1 -1  1  1
                  1 -1 -1  1]';
              
    case{'Q2'}
        degree = 2;
        dof_hat = [-1 -1 -1  0  0  0  1  1  1
                    1  0 -1 -1  0  1  1  0 -1]';
                
    case{'Q3'}
        degree = 3; 
        dof_hat = [-1  -1  -1  -1  -0.5 -0.5 -0.5 -0.5 0.5 0.5  0.5 0.5  1  1    1   1
                    1 0.5 -0.5 -1  -1   -0.5  0.5   1   1  0.5 -0.5  -1 -1 -0.5 0.5  1]';
                
                

end


dof=[];

% We compute the jacobian of the map from the reference to the physical element
% only once since we assume that all the elements have the same shape.
% If this assumption does not hold, uncomment the 3 line in the following loop,
% and modify accordingly the function map_local_physical_points.
[BJ, BJinv] = get_jacobian(region.coords_element([1:nedge],:));

for ie=1:ne
    index_el = nedge*(ie-1).*ones(1,nedge) + [1:1:nedge];
    [temp] = map_local_physical_points(region.coords_element(index_el, :), dof_hat, BJ);
    %[BJ, BJinv, temp] = map_local_physical_points(region.coords_element(index_el, :), dof_hat, Dati.type_mesh);
    dof=[dof; temp];
end

if nedge==3 % Triangles
    nln=0.5.*(degree+1).*(degree+2);
elseif nedge==4 % Rectangles
    nln=(degree+1).^2;
end

femregion=struct('fem',Dati.fem,...
    'domain',region.domain,...
    'type_mesh', region.type_mesh,...
    'nref', region.nref,...
    'h',region.h,...
    'nedges', region.nedge,...
    'nln',nln,...
    'ndof',nln*ne,...
    'ne',ne,...
    'dof',dof,...
    'nqn',Dati.nqn,...
    'degree',degree,...
    'coord',region.coord,...
    'connectivity',region.connectivity,...
    'coords_element', region.coords_element,...
    'column_elements', region.column_elements,...
    'BJ', BJ,...
    'BJinv', BJinv,...
    'boundary_edges',region.boundary_edges);

