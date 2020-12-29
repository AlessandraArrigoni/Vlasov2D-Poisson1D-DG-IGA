function [BJ, BJinv] = get_jacobian(loc_coord)
% Returns the Jacobian and its inverse of the transformation from the reference element 
% [-1,1]x[-1,1] to the current element (loc_coord).

x0=loc_coord(1,1);   % x-coordinates of vertices
x1=loc_coord(2,1);
x2=loc_coord(3,1);
x3=loc_coord(4,1);

y0=loc_coord(1,2);   % y-coordinates of vertices
y1=loc_coord(2,2);
y2=loc_coord(3,2);
y3=loc_coord(4,2);

% BJ = 0.25*[2*\delta x   0; 0  2*\delta y]
BJ(:,:) = (0.25) .* [-x0 - x1 + x2 + x3 ,  x0 - x1 - x2 + x3 ; -y0 - y1 + y2 + y3 , y0 - y1 - y2 + y3];  % Jacobian of elemental map

BJinv(:,:) = inv(BJ);
end 