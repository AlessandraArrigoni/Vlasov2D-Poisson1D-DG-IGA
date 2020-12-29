function [M] = circulant_matrix1D(femregion1D, A0, A1, trasp)
%--------------------------------------------------------------------
% PURPOSE:
%
% Assembly of the CIRCULANT matrix describing the interaction between
% neighbouring cells (following the paper Crouseilles et al. 2011).
% It is used to solve the 1D linear transport problem. 
% 
% We impose PERIODIC BOUNDARY CONDITIONS 
%
% INPUT: A0, A1 -> local matrices as defined in the paper
%        trasp -> velocity of the 1D linear transport equation (its sign will
%        determine the structure of the matrix)
%
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


fprintf('--------Begin computing CIRCULANT matrix 1D for direction %s\n',femregion1D.direction);


% Costanti per definire i vettori con indici di riga e colonna per costruire la
% matrice sparsa alla fine: pensando alle colonne.
nln = femregion1D.nln;
nmax_el = 2*nln^2; % 2*nln colonne * nln righe
nmax = femregion1D.ne*nmax_el ;
rows = zeros(nmax, 1);
cols = zeros(nmax, 1);


values_loc = reshape([A0, A1], [nmax_el,1]); % i valori sono presi per colonne
values = repmat(values_loc, [femregion1D.ne,1]);

if trasp <= 0
    for ie=1:femregion1D.ne % loop over elements
    
        % global indices of the basis functions on the current element
        index = (ie-1)*nln*ones(nln,1) + [1:nln]';
        index_neigh = ie*nln*ones(nln,1) + [1:nln]'; % indici nella cella a destra, prendendo quelli sulla prima cella se sono nell'ultima
        index_neigh(index_neigh > femregion1D.ndof) = index_neigh(index_neigh > femregion1D.ndof) - femregion1D.ndof;

        % Assemble global matrices
        rows_loc = repmat(index, [2*nln,1]); % 2 perchè per ogni cella ho sempre 2 celle con colonne non nulle
        rows((ie-1)*nmax_el + 1 : ie*nmax_el) = rows_loc;

        idx = [index ; index_neigh]'; % riga
        cols_loc = reshape(ones(nln,1)*idx, [nmax_el, 1]);
        cols((ie-1)*nmax_el + 1 : ie*nmax_el) = cols_loc;
          
    
    end

elseif trasp > 0 
    % Come sopra ma cambio le righe: ora quella che comincia con A0 non è la
    % prima ma la seconda, cioè quella descritta da index_neigh.
    for ie=1:femregion1D.ne % loop over elements
    
        % global indices of the basis functions on the current element
        index = (ie-1)*nln*ones(nln,1) + [1:nln]';
        index_neigh = ie*nln*ones(nln,1) + [1:nln]'; % indici nella cella a destra, prendendo quelli sulla prima cella se sono nell'ultima
        index_neigh(index_neigh > femregion1D.ndof) = index_neigh(index_neigh > femregion1D.ndof) - femregion1D.ndof;

        % Assemble global matrices
        rows_loc = repmat(index_neigh, [2*nln,1]); % 2 perchè per ogni cella ho sempre 2 celle con colonne non nulle
        rows((ie-1)*nmax_el + 1 : ie*nmax_el) = rows_loc;

        idx = [index ; index_neigh]'; % riga
        cols_loc = reshape(ones(nln,1)*idx, [nmax_el, 1]);
        cols((ie-1)*nmax_el + 1 : ie*nmax_el) = cols_loc;
    
    end
    
end


% Assemble global matrix
M = sparse(rows, cols, values);

end