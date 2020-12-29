%--------------------------------------------------------------------
% PURPOSE:
%
% This routine generates a structure containing all the information of the
% DG finite element space in ONE DIMENSION, i.e.,
%
% femregion.fem -> finite element space
% femregion.domain -> coordinates of the domain
% femregion.h -> mesh size
% femregion.nln -> number of local nodes
% femregion.ndof -> total number of degrees of freedom
% femregion.ne -> total number of elements
% femregion.dof -> coordinates of the degrees of freedom
% femregion.nqn -> number of 1D quadrature nodes
% femregion.degree -> polynomial degree
% femregion.coord -> coordinates of the nodes of the mesh
% femregion.coords_element -> coordinates of the elements' boundaries
%
% NB : it doesn't contain any information to link the dof 1D to the dof 2D 
% via tensor product, this is just a generic space 1D defined on the same mesh.
%
% INPUT: region -> region structure created from generate_mesh()
%        direction -> 'x' or 'y'
%--------------------------------------------------------------------

function [femregion] = create_dof1D(Data, region, direction)

fprintf('Computing dof 1D for direction %s\n',direction);

ne = sqrt(region.ne); % assuming the number of subdivisions is the same in both directions
nedge = region.nedge;

dof_hat =[];

switch Data.fem1D  % dof on the reference element [-1,1]
    case{'P0'}
        degree = 0;
        dof_hat = 0;

    case{'P1'}
        degree = 1;
        dof_hat = [-1   1]; 
        
    case{'P2'}
        degree = 2;
        dof_hat =[ -1   0   1];
               
    case{'P3'}
        degree = 3;
        dof_hat =[ -1   -0.5     0.5     1]; 

end

% Vector storing the coordinates of the elements' boundaries
coords_elem = zeros(2*ne,1); 
if direction == 'x'
    % The elements are ordered from LEFT to RIGHT
    nqn = Data.nqn(1);
    coord = region.x_points1D;
    coords_elem(1) = coord(1); 
    coords_elem(end) = coord(end);
    temp = [coord(2:end-1);
            coord(2:end-1)];
    coords_elem(2:end-1) = reshape(temp,[numel(temp),1]);
    
elseif direction == 'y'
    % The elements are ordered from BOTTOM to TOP
    nqn = Data.nqn(2);
    coord = region.y_points1D;
    coords_elem(1) = coord(1); 
    coords_elem(end) = coord(end);
    temp = [coord(2:end-1);
            coord(2:end-1)];
    coords_elem(2:end-1) = reshape(temp,[numel(temp),1]);
else
        error('Direction can only be x or y');
end
        
h = coords_elem(2) - coords_elem(1);
nln = degree + 1;
dof = zeros(nln*ne, 1);

switch Data.fem1D

    case 'P1'        
        for ie = 1:ne
            dof((ie-1)*nln*ones(nln,1) + [1:nln]') = coord(ie:ie+1); 
        end
        
    case 'P2'
        for ie = 1:ne
            dof((ie-1)*nln*ones(nln,1) + [1:nln]') = [coord(ie); coord(ie)+h/2; coord(ie+1)];
        end
        
    case 'P3'
        for ie = 1:ne
            dof((ie-1)*nln*ones(nln,1) + [1:nln]') = [coord(ie); coord(ie)+h/4; coord(ie+1)-h/4; coord(ie+1)];
        end
end



femregion = struct('fem',Data.fem1D,...
    'domain',[coords_elem(1), coords_elem(end)],...
    'direction',direction,...
    'nref', region.nref,...
    'h', h ,...
    'nln',nln,...
    'ndof',nln*ne,...
    'ne',ne,...
    'dof',dof,... % column
    'nqn',nqn,...
    'degree',degree,...
    'coord',coord,...
    'coords_element', coords_elem);

