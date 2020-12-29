function [M] = circulant_matrix1D_BC(femregion1D, A0, A1, trasp)
%--------------------------------------------------------------------
% PURPOSE:
%
% Assembly of the CIRCULANT matrix describing the interaction between
% neighbouring cells (following the paper Crouseilles et al. 2011).
% It is used to solve the 1D linear transport problem. 
% 
% We impose that the solution outside the domain is 0 (set to 0 the matrices
% referring to the neighbouring cell in this case
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
%             0     0   0    0   0 ...  A0]
%
% If trasp > 0 the global matrix has the following structure:
%           [A1    0   0    0   ... 0   0
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
rows = zeros(nmax - nln^2, 1); % Tolgo una delle matrici rispetto al caso periodico
cols = zeros(nmax - nln^2, 1);

values_loc = reshape([A0, A1], [nmax_el,1]); % i valori sono presi per colonne
values = repmat(values_loc, [femregion1D.ne - 1,1]); % Ad una cella avrò associato solo A0 o A1
values_A0 = reshape(A0, [numel(A0), 1]);
values_A1 = reshape(A1, [numel(A1), 1]);

if trasp <= 0
    
    % complete values vector
    values = [values; values_A0];

    for ie=1:femregion1D.ne-1 % loop over elements (assuming they are ordered from left to right or from bottom to top
                   
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
    
    % treat last element
    index = (femregion1D.ne-1)*nln*ones(nln,1) + [1:nln]';
    rows_loc = repmat(index, [nln, 1]);
    rows((femregion1D.ne-1)*nmax_el + 1 : end) = rows_loc;
    cols_loc = reshape(ones(nln,1)*index', [nln^2, 1]);
    cols((femregion1D.ne-1)*nmax_el + 1 : end) = cols_loc;

elseif trasp > 0 
    
    % complete values vector
    values = [values_A1; values];
        
    % treat first element: 
    index = [1:nln]';
    rows_loc = repmat(index, [nln, 1]);
    rows(1 : nln^2) = rows_loc;
    cols_loc = reshape(ones(nln,1)*index', [nln^2, 1]);
    cols(1 : nln^2) = cols_loc;
    
    % Come sopra ma cambio le righe: ora quella che comincia con A0 non è la
    % prima ma la seconda, cioè quella descritta da index_neigh.
    for ie = 2:femregion1D.ne % loop over elements (assuming they are ordered from left to right or from bottom to top
    
        % global indices of the basis functions on the current element
        index = (ie-1)*nln*ones(nln,1) + [1:nln]';
        index_neigh = ie*nln*ones(nln,1) + [1:nln]'; 
        % indici nella cella a destra, non controllo quelli fuori dai limiti perchè verranno "tirati indietro" quando costruisco le colonne

        % Assemble global matrices
        rows_loc = repmat(index, [2*nln,1]); % 2 perchè per ogni cella ho sempre 2 celle con colonne non nulle
        rows((ie-1)*nmax_el + 1 - nln^2 : ie*nmax_el - nln^2) = rows_loc; % Tolgo nln^2 perchè nella prima riga ho solo una matrice e non 2

        idx = [index ; index_neigh]' - nln; % riga : tolgo "una cella" alle colonne perchè nella matrice sono scalate indietro (non c'è A0 sulla diagonale)
        cols_loc = reshape(ones(nln,1)*idx, [nmax_el, 1]);
        cols((ie-1)*nmax_el + 1 - nln^2 : ie*nmax_el - nln^2) = cols_loc;
    
    end
    
       
end


% Assemble global matrix
M = sparse(rows, cols, values);

end