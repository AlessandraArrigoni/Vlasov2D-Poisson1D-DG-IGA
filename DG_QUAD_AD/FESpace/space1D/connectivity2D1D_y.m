function [conn] = connectivity2D1D_y(neighbours, femregion)
% PURPOSE: construct a vector with femregion.ndof entries containing the
% indices of the 2D dofs according to the folloqing order: starting from the
% bottom left of the domain and moving upwards column by column (each column is
% determined by the position of a 2D dof).
% 
% Useful to link the 2D dofs to the 1D dofs and construct a global (block diagonal)
% system with all the UNCOUPLED linear systems associated to the 1D linear
% transport problem (advection in y direction).
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
col_conn = zeros(ndof_dir, nln_dir); 

% Build LOCAL connectivity matrices to link the numbering of the dofs in each
% element to the numbering row by row that is needed. 
% We do this here since it depends only on the polynomial degree and it is 
% useless to check the condition inside the following for loop.
switch femregion.fem
    case 'Q1'
        local_numb = [2 1; 3 4]';
    case 'Q2'
        local_numb = [3 2 1; 4 5 6; 9 8 7]';
    case 'Q3'
        local_numb = [4 3 2 1; 5 6 7 8; 12 11 10 9; 13 14 15 16]';
    otherwise
        error('FEM 2D spaces can only be Q1-Q2-Q3');
end

% There are 2^nref elements per direction (the numbering is determined by the
% refining procedure used when creating the mesh)
n_el_dir = 2^femregion.nref; 
esp = 1:femregion.nref; % Start from bottom left element
cur_el = sum(2.^(2*esp-2)) + 1;
first_el_col = cur_el; % First element in the current column

for count_col = 1 : n_el_dir
    
    for count_row = 1 : n_el_dir

        % Take indices of dofs in the current element
        idx = (cur_el-1)*nln*ones(nln,1) + [1:nln]';
        
        % Sort them column by column
        idx_ord = idx(local_numb);
        col_conn((count_row-1)*nln_dir*ones(nln_dir,1) + [1:nln_dir]' , :) = idx_ord;
        
        % Move to the neighbouring element above
        cur_el = neighs(cur_el,4);
    end
    
    tot = ndof_dir*nln_dir;
    conn((count_col-1)*tot*ones(tot,1) + [1:tot]') = reshape(col_conn, [tot, 1]); % reshape takes the values columnwise
    
    % Move to the next column of elements, unless we are already considering the
    % last column on the right 
    if neighs(first_el_col, 3) == (1 + sum(2.^(2*esp-2))) % we need this condition because of the periodic BC
        break;
    else
        cur_el = neighs(first_el_col, 3);
        first_el_col = neighs(first_el_col, 3);        
    end
end


end


