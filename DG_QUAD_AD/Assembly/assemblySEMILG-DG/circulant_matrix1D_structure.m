function [rows, cols] = circulant_matrix1D_structure(femregion1D, trasp)
%--------------------------------------------------------------------
% PURPOSE:
%
% Assembly of the CIRCULANT matrix describing the interaction between
% neighbouring cells (following the paper Crouseilles et al. 2011).
% It is used to solve the 1D linear transport problem. 
% 
% We impose PERIODIC BOUNDARY CONDITIONS.
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
%            A1     0   0    0   0 ...  A0]
%
% If trasp > 0 the global matrix has the following structure:
%           [A1    0   0    0   ... 0   A0
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
rows = zeros(nmax, 1);
cols = zeros(nmax, 1);

if trasp < 0 || abs(trasp) < 1e-9
    for ie=1:femregion1D.ne % loop over elements
    
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

elseif trasp > 0 
    % As before, but we change the rows associated to the current element in the
    % global matrix: considering the first element as an example, now the 
    % "block row" starting with A0 is the second one, that is, the one linked
    % to the indices in the neighbouring element. Thus the only difference 
    % is in the construction of rows_loc.
    for ie=1:femregion1D.ne % loop over elements
    
        % Global indices of the basis functions on the current element
        index = (ie-1)*nln*ones(nln,1) + [1:nln]';
        % Global indices of the basis functions on the neighbouring cell on the
        % RIGHT (or the first cell if the current element is the last one)
        index_neigh = ie*nln*ones(nln,1) + [1:nln]'; 
        index_neigh(index_neigh > femregion1D.ndof) = index_neigh(index_neigh > femregion1D.ndof) - femregion1D.ndof;

        % Assemble vectors for global matrix 
        rows_loc = repmat(index_neigh, [2*nln,1]); % 2 since each cell interacts with another one (and ONLY one, according to our assumptions on the dt)
        rows((ie-1)*nmax_el + 1 : ie*nmax_el) = rows_loc;

        idx = [index ; index_neigh]';
        cols_loc = reshape(ones(nln,1)*idx, [nmax_el, 1]);
        cols((ie-1)*nmax_el + 1 : ie*nmax_el) = cols_loc;
    
    end
    
end


end