function [M] = global_matrixMass(femregion, rows_loc, cols_loc, values_loc)
% Builds the global block diagonal matrix collecting all the mass (or mass
% inverse) matrices for the 1D linear transport problems in one direction
% (called "local" in the following).
% The local matrices are repeated ndof1D times on the diagonal: each "block row"
% corresponds to a linear transport equation in 1D.
% 
% INPUT : femregion -> the 1D FE space associated to the set of linear transport
%                      problems
%         rows_loc, cols_loc -> vectors storing the indices of non zero rows and
%                      columns in the LOCAL matrices as given by
%                      mass_matrix1D_structure()
%         values_loc -> vector storing the values of the LOCAL mass matrix (or
%                      its inverse) as given by mass_matrix1D_structure()

dim = numel(rows_loc);
rows = zeros(femregion.ndof*dim,1);
cols = zeros(femregion.ndof*dim,1);
values = repmat(values_loc, [femregion.ndof, 1]);

for j = 1:femregion.ndof
    rows((j-1)*dim + 1 : j*dim) = (j-1)*femregion.ndof + rows_loc;
    cols((j-1)*dim + 1 : j*dim) = (j-1)*femregion.ndof + cols_loc;
end


M = sparse(rows, cols, values);
end

