function [conn] = connectivity2D1D_x(neighbours, femregion)
% PURPOSE: construct a vector with femregion.ndof entries containing the
% indices of the 2D dofs in the folloqing order: starting from the top left 
% of the rectangular domain and moving rightwards row by row (each row is
% determined by the position of a 2D dof).
%
% Useful to link the 2D dofs to the 1D dofs and construct a global (block diagonal)
% system with all the UNCOUPLED linear systems associated to the 1D linear
% transport problem (advection in x direction).
%
% INPUT : femregion -> 2D femregion
% OUTPUT : conn -> vector femregion.ndof x 1

neighs = neighbours.neigh;
conn = zeros(femregion.ndof, 1);
nln = femregion.nln;
% Only for CARTESIAN GRIDS: the number of 1D dofs (per element and total) is the
% square root of the number of 2D dofs.
nln_dir = sqrt(nln);
ndof_dir = sqrt(femregion.ndof);
row_conn = zeros(ndof_dir, nln_dir); 

% Build LOCAL connectivity matrices to link the numbering of the dofs in each
% element to the numbering row by row that is needed. 
% We do this here since it depends only on the polynomial degree and it is 
% useless to check the condition inside the following for loop.
% The local_numb matrices (before transposing them) mirror the positions of the
% dofs in the rectangular element; e.g. for Q1 the dofs are sorted as follows:
% top left, bottom left, bottom right, top right.
% We need to transpose them since the reshape (below) operates columnwise.
switch femregion.fem
    case 'Q1'
        local_numb = [1 4; 2 3]';
    case 'Q2'
        local_numb = [1 6 7; 2 5 8; 3 4 9]';
    case 'Q3'
        local_numb = [1 8 9 16; 2 7 10 15; 3 6 11 14; 4 5 12 13]';
    otherwise
        error('FEM 2D spaces can only be Q1-Q2-Q3');
end

% There are 2^nref elements per direction
n_el_dir = 2^femregion.nref;
cur_el = 1; % Start from top left element
first_el_row = 1; % First element in the current row

for count_row = 1 : n_el_dir
    
    for count_col = 1 : n_el_dir

        % Take indices of dofs in the current element
        idx = (cur_el-1)*nln*ones(nln,1) + [1:nln]';
        
        % Sort them row by row
        idx_ord = idx(local_numb);
        row_conn((count_col-1)*nln_dir*ones(nln_dir,1) + [1:nln_dir]' , :) = idx_ord;
        
        % Move to the neighbouring element on the right
        cur_el = neighs(cur_el,3);
    end
    
    tot = ndof_dir*nln_dir;
    conn((count_row-1)*tot*ones(tot,1) + [1:tot]') = reshape(row_conn, [tot, 1]); % reshape takes the values columnwise
    
    % Move to the next row of elements, unless we are already considering the
    % last row
    if neighs(first_el_row, 2) == -1 
        break;
    else
        cur_el = neighs(first_el_row, 2);
        first_el_row = neighs(first_el_row, 2);        
    end
end


end


