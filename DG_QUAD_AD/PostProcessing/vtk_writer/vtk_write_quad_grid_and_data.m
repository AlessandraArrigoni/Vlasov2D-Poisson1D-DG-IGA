function stat = vtk_write_quad_grid_and_data(filename, vtk_title, grid_X, grid_Quad, data_struct, flipped, use_lagrange_el_flag)
% modified version of vtk_write_tetrahedral_grid and data by 
% Shawn W. Walker (2016) to write 2D data on unstructured grids. Refer
% to the documentation of vtk_write_tetrahedral_grid_and_data for
% additional information

% define write commands
write_ascii = @(fid,str) [fprintf(fid,str); fprintf(fid,'\n');];
write_data  = @(fid,dat,prec) [fwrite(fid,dat,prec); fprintf(fid,'\n');];

% check the filename
[P1, N1, E1] = fileparts(filename);
if isempty(E1)
    E1 = '.vtk'; % add the extension
elseif ~strcmp(E1,'.vtk')
    disp('warning: file extension is *not* ".vtk"!');
end
filename = fullfile(P1, [N1, E1]);

% open the file in binary mode
fopen_opts = {'wb','ieee-be'};
fid = fopen(filename,fopen_opts{:});
if fid == -1
    error('Unable to write file %s: permission denied.',filename);
end

% write the initial header
write_ascii(fid,'# vtk DataFile Version 2.0');
vtk_title_reduce = vtk_title(1:min(length(vtk_title),256));
write_ascii(fid,vtk_title_reduce);
write_ascii(fid,'BINARY\n'); % give extra line-feed
% write the vertex coordinates of the mesh
write_ascii(fid,['DATASET ', 'UNSTRUCTURED_GRID']);
if ~flipped
    grid_X = grid_X'; % need to transpose for fwrite
end
GD = size(grid_X,1);
if (GD~=3)
    stat = fclose(fid);
    error('Grid vertex coordinate data is not 3-D!');
end
Num_Vtx = size(grid_X,2);
write_ascii(fid, ['POINTS ', num2str(Num_Vtx), ' float']);
write_data(fid, grid_X, 'float32');

% write the quad grid connectivity
if ~flipped
    Order = size(grid_Quad,2);
    Num_Quad   = size(grid_Quad,1);
else
    Order = size(grid_Quad,1);
    Num_Quad   = size(grid_Quad,2);
end
if (Order~=4 && ~use_lagrange_el_flag)
    stat = fclose(fid);
    error('Grid quads connectivity does not have 4 nodes per quad!');
end
Cell_Size = Num_Quad * ( Order + 1 ); % vtk needs this
write_ascii(fid,['CELLS  ', num2str(Num_Quad), '  ', num2str(Cell_Size)]);

% 1-based to 0-based indexing
min_index = min(grid_Quad(:));
if min_index==0
    disp('Quad connectivity is already 0-based!');
elseif min_index > 0
    grid_Quad = grid_Quad - 1;
else
    fclose(fid);
    error('Quad connectivity has negative indices!');
end

% modify and append extra data
if ~flipped
    DATA = uint32([(Order + 0*grid_Quad(:,1))';
                    grid_Quad']);
else
    DATA = uint32([(Order + 0*grid_Quad(1,:));
                    grid_Quad]);
end
write_data(fid,DATA, 'uint32');

% must write the cell types
% VTK has a cell type 9 for linear quads, 8 for so-called pixels ... 
write_ascii(fid,['CELL_TYPES ', num2str(Num_Quad)]);
if ( Order == 4 )
%     QUAD_LABEL = 9; 
    LABEL = 8; 
elseif(Order > 4)
    % ... and 70 for Lagrangian quads
    LABEL = 70;
else
    error('Invalid quad order!');
end
DATA = uint32(LABEL*ones(1,Num_Quad));
write_data(fid,DATA, 'uint32');

% write the POINT_DATA
if ~isempty(data_struct)
    write_ascii(fid,['POINT_DATA ', num2str(Num_Vtx)]);
    
    Num_DS = length(data_struct);
    for ii = 1:Num_DS
        if strcmpi(data_struct(ii).type,'scalar')
            write_ascii(fid,['SCALARS ', data_struct(ii).name, ' float']);
            write_ascii(fid,'LOOKUP_TABLE default');
            write_data(fid,data_struct(ii).data(:)','float32'); % scalar data easy to write
        elseif strcmpi(data_struct(ii).type,'vector')
            write_ascii(fid,['VECTORS ', data_struct(ii).name, ' float']);
            if ~flipped
                D1 = data_struct(ii).data';
            else
                D1 = data_struct(ii).data;
            end
            write_data(fid,D1,'float32');
        else
            stat = fclose(fid);
            error('Invalid data type!');
        end
    end
end

% Note: CELL_DATA is not implemented!

% end of file!
stat = fclose(fid);

end