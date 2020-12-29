function [pphys_2D] = map_local_physical_points(loc_coord, nodes_2D, BJ)
% Returns the coordinates of the 2D quadrature nodes in the current element,
% characterised by loc_coord and by the jacobian of the map BJ

x0=loc_coord(1,1);   % x-coordinates of vertices
x1=loc_coord(2,1);
x2=loc_coord(3,1);
x3=loc_coord(4,1);

y0=loc_coord(1,2);   % y-coordinates of vertices
y1=loc_coord(2,2);
y2=loc_coord(3,2);
y3=loc_coord(4,2);

% Translation vector: average of all the vertices, since the reference
% element is [-1,1]x[-1,1] and not [0,1]x[0,1]
trans = (0.25) .* [ x0 + x1 + x2 + x3 ;  y0 + y1 + y2 + y3];
        
Trasl = repmat(trans, [1, length(nodes_2D)]);
pphys_2D = (BJ*nodes_2D' + Trasl)';


end

