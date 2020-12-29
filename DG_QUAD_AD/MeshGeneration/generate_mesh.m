%--------------------------------------------------------------------
% PURPOSE:
%
% This routine generates a structure containing all the information about the
% mesh, i.e.,
%
% region.dim ->  dimension (2)
% region.type_mesh -> type of mesh
% femregion.domain -> (2x2 matrix, real) domain limits
% region.h ->  mesh size
% region.coord -> coordinates of the mesh nodes
% region.connectivity -> connectivity matrix
% region.coords_element -> coordinates of elements (counted with their multiplicity)
% region.boundary_edges -> boundary edges.
% region.column_elements -> matrix with list of elements that share the same
%                           x-interval
% Author:
% Paola Antonietti
%--------------------------------------------------------------------


function [region]=generate_mesh(Dati,step_refinement)

[g, coord, boundary_edge, connectivity]=initialize_mesh(Dati.domain, Dati.type_mesh);

switch Dati.type_mesh
    case{'TS','TU'}
        nedge=3;
    case{'CART'}
        nedge=4;
end


for i= 1:step_refinement % loop over the number of refinements
    if i==1
        switch Dati.type_mesh
            case{'TS'}
                [coord, boundary_edge, connectivity]=refine_mesh(g,coord, boundary_edge, connectivity,'regular');
                coord=jiggle_mesh(coord, boundary_edge, connectivity);
            case{'TU'}
                [coord, boundary_edge, connectivity]=refine_mesh(g,coord, boundary_edge, connectivity,'longest');
                [coord, boundary_edge, connectivity]=refine_mesh(g,coord, boundary_edge, connectivity,'longest');
                coord=jiggle_mesh(coord, boundary_edge, connectivity);
            case{'CART'}
                [coord, connectivity]= refine_quad(coord, connectivity);
        end

    else
        switch Dati.type_mesh
            case{'TS','TU'}
                [coord, boundary_edge, connectivity]=refine_mesh(g,coord, boundary_edge, connectivity,'regular');
                coord=jiggle_mesh(coord, boundary_edge, connectivity);
            case{'CART'}
                [coord, connectivity]= refine_quad(coord, connectivity);
        end
    end
end

ne=size(connectivity,2);
h=1/sqrt(ne);

coords_element=[]; % coordinates of the elements (counted with their multiplicity)
for ie=1:ne
    for k=1:nedge
        coords_element=[coords_element; coord(1,connectivity(k,ie)),coord(2,connectivity(k,ie))];
    end
end

% Find (fake) 1D elements in X and store in a matrix the corresponding
% elements in y that share the same x-interval.
% The matrix "columns" has n_intervals_x rows and n_elements_y columns: 
% associates to each 1D interval the "vertical" list of rectangles.
xpoints = unique(coord(1,:)); % it is ordered
ypoints = unique(coord(2,:));

if (strcmp(Dati.type_mesh,'CART')) 
    columns = zeros(numel(xpoints)-1, numel(xpoints)-1); % Same number of elements in x and y 
else % Original structure for triangles
    columns = zeros(numel(xpoints)-1, 2*numel(xpoints)-2); % Double number of columns (two triangles each square)
end

idx_cols = ones(numel(xpoints)-1,1); % first "free" column in each row (useful to fill it in)
for ie = 1:ne % Loop over the elements 
    xmin = min(coords_element(nedge*(ie-1)+1 : nedge*ie, 1)); 
    idx_elemx = find(xpoints == xmin);
    columns(idx_elemx, idx_cols(idx_elemx)) = ie;
    idx_cols(idx_elemx) = idx_cols(idx_elemx)+1; % new index for next iteration
end


region=struct('dim',2,...
        'type_mesh',Dati.type_mesh,...
        'nref', step_refinement,...
        'domain',Dati.domain,...
        'nedge', nedge,...
        'h',h,...
        'nvert',size(coord,2),...
        'ne',ne,...
        'coord',coord',...
        'boundary_edges',boundary_edge,...
        'connectivity',connectivity,...
        'column_elements', columns,...
        'x_points1D', xpoints,...
        'y_points1D', ypoints,...
        'coords_element',coords_element);
    
end

%--------------------------------------------------------------------
% PURPOSE:
%
% This routine initializes a tringular mesh on a given domain
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [g, coord, boundary_edge, connectivity]=initialize_mesh(domain, type_mesh)

% domain
x0=domain(1,1);
x1=domain(1,2);
y0=domain(2,1);
y1=domain(2,2);

% geometry
g=[
    2     2     2     2
    x0    x1    x1    x0
    x1    x1    x0    x0
    y1    y1    y0    y0
    y1    y0    y0    y1
    0     0     0     0
    1     1     1     1
    ];

% points: topsx - topdx - botdx - botsx
coord=[
    x0     x1     x1    x0
    y1     y1     y0    y0
    ];

% boundary edges
boundary_edge=[
    1     2     3     4
    2     3     4     1
    0     0     0     0
    1     1     1     1
    1     2     3     4
    0     0     0     0
    1     1     1     1
    ];

% connectivity
switch type_mesh
    case{'TS', 'TU'}
        connectivity =[
            1     2
            4     1
            3     3
            1     1
            ];
    case{'CART'}
        connectivity =[
            1
            4
            3
            2
            1
            ];
end

end