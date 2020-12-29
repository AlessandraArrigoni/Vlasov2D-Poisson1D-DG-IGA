function [dphiq, coords] = write_solDG_to_vtk (filename, Data, femregion, u_h, npoints, plot_perturb, dphiq, coords)
% Write to .vtk the solution on each SQUARE element by evaluating the basis functions on
% a given set of points (not necessarily the quadrature nodes).
% We need this because we want to use high order polynomials and the matlab
% function patch connects the values in a linear way.

% INPUT : filename = the .vtk output file
%         u_h = the set of coefficients with respect to the DG basis
%         npoints = number of points [nx,ny] where we want to evaluate the
%                   solution (the first and last one are on the boundaries of the elem)
%         plot_perturb = boolean to add the perturbation evolution to the output
%                   file or not (equilibrium is Data.f_eq)
%         dphiq = {optional} shape functions already evaluated at the given set
%                   of points
%         coords = {optional} coordinates of all the points in the whole mesh
%                   (in case we already computed them for a previous time) 

nln=femregion.nln;
ne=femregion.ne;

if nargin < 7
    % Points on the reference element [-1,1]^2
    xref = linspace(-1,1,npoints(1));
    yref = linspace(-1,1,npoints(2));
    points2D = []; % They have the same structure as node_2D from quadrature()
    for i = 1:npoints(1)
        for j = 1:npoints(2)
            points2D = [points2D; xref(i), yref(j)];
        end
    end

    % scalar shape functions
    [shape_basis]= basis_lagrange(Data.fem);

    % evaluation of shape functions
    [dphiq, ~, ~]= evalshape(femregion, shape_basis, points2D, {xref, yref}, femregion.nln);
    dphiq = reshape(dphiq,[length(points2D),femregion.nln]); % nqn2D x nln
    
    coords = []; 
end

% Solution in all the points of all elements
sol = []; 

for ie=1:ne % loop over elements
    
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    index_element=femregion.nedges*(ie-1).*ones(femregion.nedges,1) + [1:1:femregion.nedges]';

    local_uh = u_h(index);
    local_sol = dphiq * local_uh; % Column vector n_points2D long 
    
    sol = cat(3,sol,local_sol);
   
    if nargin < 7
        % translate the points on the current element and create the grid 
        coords_elem=femregion.coords_element(index_element, :);
        [pphys_2D] = map_local_physical_points(coords_elem, points2D, femregion.BJ); %  Jacobian and physical coordinates at the quadrature points

        
        pphys_2D = reshape(pphys_2D', [2, npoints(2), npoints(1)]); % 2 is the domain dimension (2D)
        coords = cat(4, coords, pphys_2D);
    end
    
end

data_title = 'solution plot';

x = reshape(coords(1,:,:,:),[],1); 
y = reshape(coords(2,:,:,:),[],1); 
z = zeros(size(x)); % Fake z coordinate, needed for vtk file

data_struct(1).type = 'scalar';
u = reshape(sol,[],1);
        
data_struct(1).name = 'solution';
data_struct(1).data = u';

if plot_perturb
    data_struct(2).type = 'scalar';
    data_struct(2).name = 'perturbation';
    df = u - Data.f_eq(x,y);
    data_struct(2).data = df';
end

% Construct connectivity for each element: for each cell the order of the
% vertices required by vtk is bottom left - bottom right - top left - top right
nel_dir = npoints - 1; % number of "cells" per direction in each element
cells_per_el = prod(nel_dir);
connectivity_el = zeros(4, cells_per_el);
points_per_el = prod(npoints);
m = npoints(2);

for i = 1:npoints(1)-1
    index = i;
    
    for j = 1: m-1
        
        connectivity_el(1,index) = j + m*(i-1);
        connectivity_el(2,index) = j + m*i;
        connectivity_el(3,index) = j+1 + m*(i-1);
        connectivity_el(4,index) = j+1 + m*i;
        
        index = index + nel_dir(1);
    end
end

connectivity = zeros(4,cells_per_el*ne);
idx = 1; shifting = 0;
for i = 1:ne
    connectivity(:, idx :(idx + cells_per_el - 1)) = connectivity_el + shifting;
    idx = idx + cells_per_el; 
    shifting = shifting + points_per_el;
end

data_points = [x y z]';

vtk_write_quad_grid_and_data(filename, data_title, data_points, connectivity, data_struct, true);


end