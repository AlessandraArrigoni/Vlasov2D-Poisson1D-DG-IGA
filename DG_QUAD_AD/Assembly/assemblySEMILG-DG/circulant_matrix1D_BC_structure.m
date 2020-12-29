function [rows, cols] = circulant_matrix1D_BC_structure(femregion1D, trasp)
%--------------------------------------------------------------------
% PURPOSE:
%
% Assembly of the CIRCULANT matrix describing the interaction between
% neighbouring cells (following the paper Crouseilles et al. 2011).
% It is used to solve the 1D linear transport problem. 
% 
% We impose "COMPACT SUPPORT" BOUNDARY CONDITIONS: the solution outside the
% domain is 0, that is, the local matrix referring to the neighbouring cell
% is set to 0.
% Only the rows and columns structure to build the sparse matrix is returned
% since it is always the same (according to the SIGN of the transport term)
%
% INPUT: trasp -> velocity of the 1D linear transport equation (its sign will
%        determine the structure of the matrix)
% 
% If trasp <= 0 the global matrix has the following structure (with A1 = 0 if trasp = 0):
%           [A0    A1   0    0   0 ...  0
%             0    A0   A1   0   0 ...  0
%           ...                    ...
%             0     0   0    0   0 ...  A0]
%
% If trasp > 0 the global matrix has the following structure:
%           [A1    0   0    0   ... 0   0
%            A0    A1  0    0   ... 0   0
%           ...                 ... 
%            0     0   0    0   ... A0  A1]
%--------------------------------------------------------------------------


% fprintf('--------Begin computing STRUCTURE CIRCULANT matrix 1D for direction %s\n',femregion1D.direction);

% Useful constants to build the vectors with the rows and columns' indices
% needed to assemble the sparse matrix at the end of the function.
nln = femregion1D.nln;
nmax_el = 2*nln^2; % 2*nln columns * nln rows
nmax = femregion1D.ne*nmax_el ;
% We need one matrix less compared to the periodic BC case
rows = zeros(nmax - nln^2, 1); 
cols = zeros(nmax - nln^2, 1);

if trasp < 0 || abs(trasp) < 1e-9

    for ie=1:femregion1D.ne-1 % loop over elements (until the last but one)
        
        % Global indices of the basis functions on the current element
        index = (ie-1)*nln*ones(nln,1) + [1:nln]';
        % Global indices of the basis functions on the neighbouring cell on the
        % RIGHT (or the first cell if the current element is the last one)
        index_neigh = ie*nln*ones(nln,1) + [1:nln]'; 
        index_neigh(index_neigh > femregion1D.ndof) = index_neigh(index_neigh > femregion1D.ndof) - femregion1D.ndof;
        
        % Assemble vectors for global matrix 
        rows_loc = repmat(index, [2*nln,1]); % 2 since each cell interacts with another one (and ONLY one, according to our assumptions on the dt)
        rows((ie-1)*nmax_el + 1 : ie*nmax_el) = rows_loc;

        idx = [index ; index_neigh]';
        cols_loc = reshape(ones(nln,1)*idx, [nmax_el, 1]);
        cols((ie-1)*nmax_el + 1 : ie*nmax_el) = cols_loc;
              
    end
    
    % Treat last element (no columns associated to the "neighbouring" element)
    index = (femregion1D.ne-1)*nln*ones(nln,1) + [1:nln]';
    rows_loc = repmat(index, [nln, 1]);
    rows((femregion1D.ne-1)*nmax_el + 1 : end) = rows_loc;
    cols_loc = reshape(ones(nln,1)*index', [nln^2, 1]);
    cols((femregion1D.ne-1)*nmax_el + 1 : end) = cols_loc;

elseif trasp > 0 
    
    % Treat first element (no columns associated to the "neighbouring" element)
    index = [1:nln]';
    rows_loc = repmat(index, [nln, 1]);
    rows(1 : nln^2) = rows_loc;
    cols_loc = reshape(ones(nln,1)*index', [nln^2, 1]);
    cols(1 : nln^2) = cols_loc;
    
    % As before, but we change the columns associated to the current element in the
    % global matrix: now the matrix A0 (on the current element) is no longer on
    % the diagonal, but it is on the position linked to the "previous" element.
    % Thus we change the "cols_loc" vector by subtracting nln (number of local
    % basis functions) at each time.
    for ie = 2:femregion1D.ne % loop over elements (from the second one)
        
        % Global indices of the basis functions on the current element
        index = (ie-1)*nln*ones(nln,1) + [1:nln]';
        index_neigh = ie*nln*ones(nln,1) + [1:nln]'; 
        
        % Assemble vectors for global matrix 
        rows_loc = repmat(index, [2*nln,1]); % 2 since each cell interacts with another one (and ONLY one, according to our assumptions on the dt)
        % We must subtract nln^2 since in the first "block row" there is only
        % one matrix (not two), so the number of elements in the "rows" vector
        % is reduced as well.
        rows((ie-1)*nmax_el + 1 - nln^2 : ie*nmax_el - nln^2) = rows_loc; 
        
        idx = [index ; index_neigh]' - nln; 
        cols_loc = reshape(ones(nln,1)*idx, [nmax_el, 1]);
        cols((ie-1)*nmax_el + 1 - nln^2 : ie*nmax_el - nln^2) = cols_loc;
    
    end
    
end


end