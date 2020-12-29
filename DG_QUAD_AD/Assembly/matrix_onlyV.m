%--------------------------------------------------------------------
% PURPOSE:
%
% Assembly of the matrix corresponding to the transport operator with jump
% stabilization (upwinding technique).
% Only the terms containing the FIRST COMPONENT of the transport vector beta 
% are assembled here (volume term and integrals on the vertical edges 1, 3) 
% since they do not depend on time.
% 
% ADV = sparse(femregion.ndof,femregion.ndof); \sum(- \int_{element} beta * u * grad(v) dx)
% ED = sparse(femregion.ndof,femregion.ndof);  \sum(\int_{element_boudary} {beta u}_beta * [v] ds)
%
% The rhs is not computed here since it depends on time.
%
% The number of quadrature nodes can be different in x and y direction
%
%--------------------------------------------------------------------


function [L]= matrix_onlyV(femregion,neighbour,Data,basis)


% Extract the fields from the data structures into local variables
data_names = fieldnames (basis);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= basis.(data_names{iopt});']);
end
data_names = fieldnames (Data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= Data.(data_names{iopt});']);
end

% Only the 1D nodes in the y direction are needed here.
w_1D_y = w_1D{2};
nqn_y = length(w_1D_y);

% Get Jacobian: assuming all the elements have the same shape
nln = femregion.nln;
BJ = femregion.BJ;
Jdet = abs(det(BJ));
BJinv = femregion.BJinv;

% Useful constants to build the vectors with the rows and columns' indices
% needed to assemble the sparse matrix at the end of the function.
% nmax_el: in each element, each trial function interacts at most with all the
% test functions in the same element (nln) and with all the test functions in
% one of the neighbouring element (nln), according to the upwind principle and
% to the assumption that the transport coefficient cannot change sign in the
% same element.
nmax_el = 2*nln^2; % 2*nln columns * nln rows
nmax = femregion.ne*nmax_el ;
rows = zeros(nmax, 1);
cols = zeros(nmax, 1);
values = zeros(nmax, 1);
last_idx = 0; % Index of the current last element in the global vectors.


for ie=1:femregion.ne % loop over elements
       
    % global indices of the basis functions on the current element
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';

    index_element=femregion.nedges*(ie-1).*ones(femregion.nedges,1) + [1:1:femregion.nedges]';

    neigh_ie=neighbour.neigh(ie,:);
    neighedges_ie=neighbour.neighedges(ie,:);

    coords_elem=femregion.coords_element(index_element, :);

    [pphys_2D] = map_local_physical_points(coords_elem, nodes_2D, BJ);
    [~, ~, pphys_1D] = get_jacobian_physical_points_faces(coords_elem, nodes_1D, type_mesh);
    % pphys_1D is a cell array with one cell per face
               
    % We assume that the mesh is structured and all the elements have the same
    % shape -> we can compute the normals and the meshsize only once
    if ie == 1 
        [normals,meshsize] = get_normals_meshsize_faces(coords_elem);
    end
    
    ADV_loc = zeros(femregion.nln, femregion.nln);
 
    dx = w_2D*Jdet;
    Y = pphys_2D(:,2); % Discrete values of the transport coefficient (according to Vlasov equation)
    
    for k=1:length(w_2D) % loop over 2D quadrature nodes
        ADV_loc = ADV_loc - ([Y(k),0]*(BJinv'* gradphi(:,:,k)))'*dphiq(k,:)*dx(k);        
    end 

    ED_loc = zeros(femregion.nln, femregion.nln);

    j = [];
    EDN = []; % Edge terms on the neighbouring side. 

    for iedg = [1,3] % Only vertical sides: 1 = left, 3 = right
        EDNtemp = zeros(femregion.nln, femregion.nln);
        neigedge = neighedges_ie(iedg);    % index of neighbour edge
        ds = meshsize(iedg)*w_1D_y; 
        Y = pphys_1D{iedg}(:,2); % Discrete values of the transport coefficient (according to Vlasov equation)
        
        for k=1:nqn_y   % loop over 1D quadrature nodes            
            
            kk = nqn_y + 1 - k;  % index of neighbouring quad point            
            traspE = [Y(k),0];
           
            ED_loc = ED_loc ...                          
                 + 0.5.*(traspE*normals(:,iedg) + abs(traspE*normals(:,iedg))).*B_edge{iedg}(:,k).*B_edge{iedg}(:,k)'.*ds(k);
                       
            if neigh_ie(iedg) ~= -1 %internal faces
                
                %  \int_{E_h} {beta u}_beta * [v] ds
                EDNtemp  = EDNtemp ...
                     + 0.5.*(traspE*normals(:,iedg) - abs(traspE*normals(:,iedg))).* B_edge{iedg}(:,k).*B_edge{neigedge}(:,kk)'.*ds(k);
            end
            
                       
        end
        if neigh_ie(iedg) ~= -1 %internal faces
            j = [j; (neigh_ie(iedg)-1)*nln*ones(nln,1) + [1:nln]'];
            EDN = [EDN, EDNtemp];
        end
    end
    
   
    % Assemble vectors for global matrix 
    size_loc = nln*(nln+numel(j)); % number of non zero entries on the current element
    values_loc1 = [ED_loc + ADV_loc, EDN];
    values_loc = reshape(values_loc1, [size_loc,1]);
    values(last_idx + 1 : last_idx + size_loc) = values_loc;
    
    rows_loc = repmat(index, [nln+numel(j),1]);
    rows(last_idx + 1 : last_idx + size_loc) = rows_loc;
    
    indices = [index; j]';
    cols_loc = reshape(ones(nln,1)*indices, [nln*(nln+numel(j)),1]);
    cols(last_idx + 1 : last_idx + size_loc) = cols_loc;
 
    last_idx = last_idx + size_loc;
    
end

% Assemble global matrix
L = sparse(rows(1:last_idx), cols(1:last_idx), values(1:last_idx));

end

