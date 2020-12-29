%% TEST SOLID BODY ROTATION with SEMILAGRANGIAN DG METHOD
% Script to test the implementation of the SEMILG-DG method on the linear 
% transport equation 2D with II order splitting.
% We apply periodic in x and "compact support" (i.e. the solution
% is 0 outside the domain) in y boundary conditions.
%
% We compute the errors in L2 norm with respect to the exact solution.
% NB: the timestep must be sufficiently small to ensure that the
% characteristics' starting point falls in the neighbouring cell at most for
% every possible value of the two components of the transport vector. The 
% timestep for advection in v is dt/2 so we can choose dt_scalings larger than
% expected.
%
% To understand the meaning of "CIRCULANT MATRIX" refer to the paper by
% Crouseilles et al. (2011)
% 
% The problem we consider is: df/dt + u(y)* df/dx + u(x)* df/dy = 0
% The exact solution features a smooth profile rotating around (x,y) = (0.4,0.5)
% as in the paper by Rossmanith & Seal (2011)

data = struct('domain',      [0,1; 0,1],...
              'initial_f',   @(r,t) cos(5*pi*r/3).^6.*(r<=0.3),... % t is useless in this problem, but needed by the postprocessing
              'test',        'solid_body',...
              'x0',          0.4,...
              'y0',          0.5,...
              'fem',         'Q2',...
              'fem1D',       'P2',...
              'BC',          'Period',... 
              'nqn',         [3,3],...  % Number of quadrature nodes in (x,v)
              'time',         0,...
              'Tend',         1,...
              'type_mesh',    'CART');
          
nref = [3,4,5,6,7];     % Refinement levels: 2^nref cells per direction
scalings = [4,8,12,16]; % It must be dt < dx/trasp 
hh = [];

[nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(data.nqn);
% Translate nodes on the reference element (0,1). We can do it only once if we
% assume that we set the same number of nodes in both directions.
nodes_ref = (nodes_1D{1} + 1)*0.5; 

% REFINEMENTS LOOP
count = 0;
for ref = nref % Loop over refinements
    fprintf(['\nSTART WITH REF ' num2str(ref) '\n']);
    count = count + 1;
    errL2 = [];
    
    % MESH AND SPACE FOR DG METHOD
    region = generate_mesh(data, ref);    % Region 2D
    femregion = create_dof(data, region); % FE space 2D
    femregionX = create_dof1D(data, region, 'x'); % FE space 1D in x
    femregionY = create_dof1D(data, region, 'y'); % FE space 1D in y
    if ref < 4
        neighbour = neighbours_Period(femregion); 
    else % load saved neighbour structure to save computational time
        load(['neighbour_ref' num2str(ref) '.mat']);
    end
    % Build connectivity vectors to link 2D dofs to 1D dofs
    connX = connectivity2D1D_x(neighbour, femregion);
    connY = connectivity2D1D_y(neighbour, femregion);
    
    % Evaluate basis 2D on the REFERENCE ELEMENT (-1,1)x(-1,1)
    shape2D = basis_lagrange(femregion.fem); % scalar shape functions
    [dphiq, Grad, B_edge] = evalshape(femregion, shape2D, nodes_2D, nodes_1D, femregion.nln);
    grad2D = permute(Grad,[2,3,1]); % 2 x nln x nqn2D
    phi2D = reshape(dphiq,[length(nodes_2D), femregion.nln]); % nqn2D x nln    
    basis2D = struct('nodes_1D', {nodes_1D}, 'w_1D', {w_1D}, 'nodes_2D', nodes_2D, 'w_2D', w_2D, ...
            'dphiq', phi2D, 'gradphi', grad2D, 'B_edge', {B_edge});
        
    % Evaluate basis 1D on the REFERENCE ELEMENT (0,1)
    % We can do it only once for X and Y since we have the same number of 
    % quadrature nodes.
    shape1D = basis_lagrange1D(data.fem1D);    
    [phi1D, grad1D] = evalshape1D(shape1D, nodes_ref, femregionX.nln); 
    basis1D = struct('nodes_1D', nodes_1D{1}, 'w_1D', w_1D{1}, ...
                'dphiq', phi1D, 'gradphi', grad1D);

    % MASS MATRIX 1D (inverse): local and global for the 2 directions; the
    % values and the structure is the same for all the uncoupled 1D linear
    % systems in each direction, but they differ in the two directions.
    [rowsX, colsX, valuesX, valuesinvX] = mass_matrix1D_structure(femregionX, basis1D);
    MinvX = global_matrixMass(femregionX, rowsX, colsX, valuesinvX);
    
    [rowsY, colsY, valuesY, valuesinvY] = mass_matrix1D_structure(femregionY, basis1D);
    MinvY = global_matrixMass(femregionY, rowsY, colsY, valuesinvY);

    % Local CIRCULANT MATRIX structures: we precompute them since they only
    % depend on the SIGN of the transport coefficient; thus we consider both
    % cases (negative/0 or positive) through the last input parameter +1/-1
    % Periodic BC
    [rows_loc_plus, cols_loc_plus] = circulant_matrix1D_structure(femregionX, +1);
    [rows_loc_minus, cols_loc_minus] = circulant_matrix1D_structure(femregionX, -1);
    structures_circX = struct('rows_loc_plus', rows_loc_plus, 'cols_loc_plus',cols_loc_plus,...
                            'rows_loc_minus', rows_loc_minus, 'cols_loc_minus', cols_loc_minus);
    % "Compact support" BC
    [rows_loc_plus, cols_loc_plus] = circulant_matrix1D_BC_structure(femregionY, +1);
    [rows_loc_minus, cols_loc_minus] = circulant_matrix1D_BC_structure(femregionY, -1);
    structures_circY = struct('rows_loc_plus', rows_loc_plus, 'cols_loc_plus',cols_loc_plus,...
                            'rows_loc_minus', rows_loc_minus, 'cols_loc_minus', cols_loc_minus);
    
    % BUILD structs storing the TRANSPORT VALUES
    % unique_coord_dof is a cell array with the indices of the 1D dofs
    % associated to the corresponding value of the transport coefficient (stored
    % in "values"). The numbering of the 1D dofs is needed to identify the local
    % systems (one per value) into the global matrix collecting all of them.
    % We are not interested in the "real" numbering of the dofs, but in the
    % order of the linear systems: we choose to start from the top, so we must
    % flip the coordinates stored in femregionY
    dof1DY = flip(unique(femregionY.dof));
    temp = flip(femregionY.dof);
    for i = 1:length(dof1DY)
        unique_coord_dofY{i} = find( temp == dof1DY(i));
    end
    uY = pi*(2*dof1DY - 1);
    Vstruct = struct('values', uY, 'unique_coord_dof', {unique_coord_dofY});
    
    % Helper structure to deal with E in building the matrix for the transport
    % in x direction.
    % In the Vlasov-Poisson case we cannot precompute the values uX since they 
    % depend on time and on the previous solutions, so, to keep the same code
    % structure, we choose not to compute them even in this case, although the
    % values are constant.
    dof1DX = unique(femregionX.dof);
    for i = 1:length(dof1DX)
        unique_coord_dofX{i} = find( femregionX.dof == dof1DX(i));
    end
    
    for dt_scaling = scalings   % Loop over dt_scalings
               
        % We take the smallest dx to find the dt satisfying our assumption on
        % the starting point for the characteristics.
        dx = min(femregionX.h, femregionY.h); 
        dt = dx/dt_scaling;
        
        % Initialize the system
        rvett = sqrt((femregion.dof(:,1)-data.x0).^2+(femregion.dof(:,2)-data.y0).^2);
        f_old = data.initial_f(rvett,0); 
        f_quarter = zeros(size(f_old));
        f_half = zeros(size(f_old));
        f_new = zeros(size(f_old));
        
        % BUILD GLOBAL MATRIX for the transport in X direction, collecting all
        % the systems for the 1D transport problems associated to the discrete
        % values of the transport coefficient.
        Av = global_matrixV(femregionX, femregionY, shape1D, nodes_ref, w_1D{1}, 0.5*dt, structures_circX, Vstruct); 
    
        % Time loop
        for t = 0 : dt : data.Tend-dt

            % ADVECTION IN X (HALF TIME STEP) 
            temp = MinvX * Av * f_old(connX); % f_old(connX) changes the order of the dofs, from 2D to the chosen 1D ordering 
            f_quarter(connX) = temp; % temp contains the values of the dofs ordered in 1D -> we want them in the 2D order

            % BUILD structs storing the TRANSPORT VALUES
            uX = pi*(1-2*dof1DX);
            Estruct = struct('values', uX, 'unique_coord_dof', {unique_coord_dofX});

            % ADVECTION IN V (FULL TIME STEP)
            Ae = global_matrixE(femregionX, femregionY, shape1D, nodes_ref, w_1D{1}, dt, structures_circY, Estruct); 
            temp = MinvY * Ae * f_quarter(connY);
            f_half(connY) = temp;

            % ADVECTION IN X (HALF TIME STEP) 
            temp = MinvX * Av * f_half(connX);
            f_new(connX) = temp;

            f_old = f_new;

        end
    
        % COMPUTE ERRORS AT FINAL TIME
        data.time = t+dt;
        errors = compute_errors(data, femregion, f_new);
        errL2 = [errL2 errors.E_L2];
    
    end % Loop over dt scalings
    
    errors_ref{count} = struct('ref', ref, 'errL2', errL2);
    hh = [hh, femregion.h];
    fprintf(['FINISH REF ' num2str(ref) '\n'])
    
end % Loop over refinements

fprintf('!!! END TEST !!!\n')

%% Convergence analysis
% (it should be modified if more than one dt_scaling is tested)

figure()
loglog( hh, errL2, 'o-', hh, hh.^(femregionX.degree -1 ),'--')
legend('L2 norm',['h^' num2str(femregionX.degree-1 )]);
title(['Errors with polynomials ' femregionX.fem])


pL2 = log(errL2(2:end)./errL2(1:end-1))./log(hh(2:end)./hh(1:end-1))

%% RESULTS ANALYSIS WITH DIFFERENT CFL CONDITIONS
errors = zeros(length(scalings), length(nref));
hh2D = hh;

for j = 1:length(nref)
    for  i = 1:length(scalings)
        errors(i,j) = errors_ref{j}.errL2(i);
    end
end

% Plot errors
colors = hsv(length(scalings));
markers = ['o','*','d','+']; % Knowing that we used 4 different scalings
figure()
for i = 1:length(scalings)
    loglog(hh2D, errors(i,:),'-','Marker',markers(i),'Color',colors(i,:),'LineWidth',1)
    hold on
    leg{i} = ['CFL = \pi/' num2str(scalings(i))];
    
end

if strcmp(data.fem,'Q2')
    % CASE Q2
    loglog(hh2D, 0.3*hh2D.^(femregion.degree),'k--', hh2D, hh2D.^(femregion.degree + 1),'k-.','LineWidth',1); 
    leg{i+1} = ['h^' num2str(femregion.degree)];
    leg{i+2} = ['h^' num2str(femregion.degree+1)];
    legend(leg,'Location','SE','FontSize',14); xlabel('h','FontSize',14); ylabel(['Error L^2 at T = ' num2str(data.Tend)],'FontSize',14);
    title(['Errors with polynomials ' data.fem],'FontSize',18) 
    
elseif strcmp(data.fem,'Q3')
    % CASE Q3
    loglog(hh2D, 0.5*hh2D.^(femregion.degree-1),'k--', hh2D, 0.5*hh2D.^(femregion.degree),'k-.','LineWidth',1); 
    leg{i+1} = ['h^' num2str(femregion.degree-1)];
    leg{i+2} = ['h^' num2str(femregion.degree)];
    legend(leg,'Location','SE','FontSize',14); xlabel('h','FontSize',14); ylabel(['Error L^2 at T = ' num2str(data.Tend)],'FontSize',14);
    title(['Errors with polynomials ' data.fem],'FontSize',18) 
end

