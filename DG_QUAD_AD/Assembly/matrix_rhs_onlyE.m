%--------------------------------------------------------------------
% PURPOSE:
%
% Assembly of the source term and of the matrix corresponding to the transport
% operator with jump stabilization (upwinding technique).
% Only the terms containing the SECOND COMPONENT of the transport vector beta 
% are assembled here (volume term and integrals on the vertical edges 2, 4) 
% since they depend on time.
% 
% ADV = sparse(femregion.ndof,femregion.ndof); \sum(- \int_{element} beta * u * grad(v) dx)
% ED = sparse(femregion.ndof,femregion.ndof);  \sum(\int_{element_boudary} {beta u}_beta * [v] ds)
% f = vector(femregion.ndof); \sum(\int_{element} f * v dx dv)
%
% The INPUT beta_vett is a matrix (n_quad_nodes_per_el X n_el in x direction)
% with the values of the electric field associated to each 1D node in x. 
% 
% The UPWIND FLUX is determined according to the sign of the CELL AVERAGE of E: 
% since it is time-dependent it may change sign inside a cell iteration by
% iteration.
% 
% The variable "time" is used just for the rhs since the
% transport does not depend explicitly on t.
%
%--------------------------------------------------------------------


function [L, f]= matrix_rhs_onlyE(femregion,neighbour,Data,basis,beta_vett)


% Extract the fields from the data structures into local variables
data_names = fieldnames (basis);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= basis.(data_names{iopt});']);
end
data_names = fieldnames (Data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= Data.(data_names{iopt});']);
end

% Only the 1D nodes in the x direction are needed here.
w_1D_x = w_1D{1};
nqn_1D_x = length(w_1D_x);
nqn_1D_y = length(w_1D{2});
    
% Get Jacobian: assuming all the elements have the same shape
nln = femregion.nln;
BJ = femregion.BJ;
Jdet = abs(det(BJ));
BJinv = femregion.BJinv;

% Compute "CELL AVERAGES" for E (beta_vett): \int_{cell} E_h dx
% Since we are only interested in its sign we forget the scaling with meshsize.
sign_avg = sign(w_1D_x * beta_vett); % Vector with the signs of the "average" for each cell in x. 

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

f = zeros(femregion.ndof, 1);

for ie=1:femregion.ne % loop over elements
    
    % Find index of 1D element that shares the same x-interval as the current 2D
    % element ie.
    % femregion.column_elements has n_elements1D rows storing the indices of the
    % corresponding 2D elements
    [idx_1D, ~] = find(femregion.column_elements == ie);
    
    % Global indices of the basis functions on the current element
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';

    index_element=femregion.nedges*(ie-1).*ones(femregion.nedges,1) + [1:1:femregion.nedges]';

    neigh_ie=neighbour.neigh(ie,:);
    neighedges_ie=neighbour.neighedges(ie,:);

    coords_elem=femregion.coords_element(index_element, :);

    [pphys_2D] = map_local_physical_points(coords_elem, nodes_2D, BJ);
        
    % We assume that the mesh is structured and all the elements have the same
    % shape -> we can compute the normals and the meshsize only once
    if ie == 1
        [~,meshsize] = get_normals_meshsize_faces(coords_elem);
    end
    
    ADV_loc = zeros(femregion.nln, femregion.nln);
    f_loc = zeros(femregion.nln,1);
    
    % Evaluation of the source function on the quadrature nodes
    F = source(pphys_2D(:,1), pphys_2D(:,2),time);
    
    dx = w_2D*Jdet;
    
    for k=1:length(w_2D) % loop over 2D quadrature nodes

        % Select row corresponding to the 1D node that has the same x coordinate
        % as the current 2D node.
        row = floor(k/nqn_1D_y) + 1*(mod(k,nqn_1D_y)>0);
        trasp = [0, beta_vett(row, idx_1D)];
        
        f_loc = f_loc + F(k)*dphiq(k,:)'*dx(k);
        
        ADV_loc = ADV_loc - (trasp*(BJinv'* gradphi(:,:,k)))'*dphiq(k,:)*dx(k);
        
    end % End loop over 2D nodes        
   
    ED_loc = zeros(femregion.nln, femregion.nln);
    
    j = [];
    EDN = [];  % Edge terms on the neighbouring side. 
    EDNtemp = zeros(femregion.nln, femregion.nln);
    
    % Only horizontal sides: the code is repeated since we must assemble
    % different terms according to the side.
    % SIDE 2 (bottom)
    iedg = 2;
    neigedge = neighedges_ie(iedg);    % index of neighbour edge
    betaE = beta_vett(:,idx_1D);       % nodes sorted from left to right 
    ds = meshsize(iedg)*w_1D_x; 
    
    if sign_avg(idx_1D) >= 0 % check on the average sign
        % Assemble neighbour or boundary term (the normal is [0,-1])
        if neigh_ie(iedg) ~= -1 % internal faces
            j = [j; (neigh_ie(iedg)-1)*nln*ones(nln,1) + [1:nln]'];
            
            for k=1:nqn_1D_x    % loop over 1D quadrature nodes
            
                kk = nqn_1D_x + 1 - k;  % index of neighbouring quad point
          
                EDNtemp = EDNtemp - ...
                    betaE(k).* B_edge{iedg}(:,k).*B_edge{neigedge}(:,kk)'.*ds(k);
            end
            
            EDN = [EDN, EDNtemp];
        end
        
    elseif sign_avg(idx_1D) < 0 
        % Assemble internal term on the current element 
        if neigh_ie(iedg) ~= -1 % internal faces
            for k=1:nqn_1D_x    % loop over 1D quadrature nodes
                        
                ED_loc = ED_loc - betaE(k).*B_edge{iedg}(:,k).*B_edge{iedg}(:,k)'.*ds(k);
                
            end
        end
    end
    
    % SIDE 4 (TOP)
    iedg = 4;
    neigedge = neighedges_ie(iedg);           % index of neighbour edge
    betaE = swap_vector(beta_vett(:,idx_1D)); % nodes sorted from rigth to left 
    ds = meshsize(iedg)*w_1D_x; 
    
    EDNtemp = zeros(femregion.nln, femregion.nln);
    
    if sign_avg(idx_1D) < 0 % check on the average sign
        % Assemble neighbour or boundary term (the normal is [0,1])
        if neigh_ie(iedg) ~= -1 % internal faces
            j = [j; (neigh_ie(iedg)-1)*nln*ones(nln,1) + [1:nln]'];
            
            for k=1:nqn_1D_x    % loop over 1D quadrature nodes
        
                kk = nqn_1D_x + 1 - k;  % index of neighbouring quad point
                
                EDNtemp = EDNtemp ...
                    + betaE(k).* B_edge{iedg}(:,k).*B_edge{neigedge}(:,kk)'.*ds(k);

            end
            EDN = [EDN, EDNtemp];
        end
              
    elseif sign_avg(idx_1D) >= 0 
        % Assemble internal term on the current element 
        if neigh_ie(iedg) ~= -1 % internal faces
            for k=1:nqn_1D_x    % loop over 1D quadrature nodes
                
                ED_loc = ED_loc + betaE(k).*B_edge{iedg}(:,k).*B_edge{iedg}(:,k)'.*ds(k);

           end
        end
    end
    
    % Assemble rhs
    f(index) = f(index) + f_loc; 
   
    % Assemble vectors for global matrix 
    size_loc = nln*(nln+numel(j)); % number of non zero entries on the current element
    values_loc1 = [ED_loc + ADV_loc, EDN];
    values_loc = reshape(values_loc1, [size_loc,1]);
    values(last_idx + 1 : last_idx + size_loc) = values_loc;
   
    rows_loc = repmat(index, [nln+numel(j),1]);
    rows(last_idx + 1 : last_idx + size_loc) = rows_loc;

    
    indices = [index; j]';
    cols_loc = reshape(ones(nln,1)*indices, [size_loc,1]);
    cols(last_idx + 1 : last_idx + size_loc) = cols_loc;
    
    last_idx = last_idx + size_loc;
    
end


% Assemble global matrix
L = sparse(rows(1:last_idx), cols(1:last_idx), values(1:last_idx));

end

