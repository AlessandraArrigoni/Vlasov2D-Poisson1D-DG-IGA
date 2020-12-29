%% CONVERGENCE TEST with SEMILAGRANGIAN DG METHOD
% Script to test the implementation of the SEMILG-DG method on the Vlasov-Poisson
% system in 2D with II order splitting.
% We apply periodic in x and "compact support" (i.e. the solution
% is 0 outside the domain) in y boundary conditions.
%
% We compute the errors in L2 norm with respect to the exact solution.
% NB: the timestep must be sufficiently small to ensure that the
% characteristics' starting point falls in the neighbouring cell (at most) for
% every possible value of the two components of the transport vector. The 
% timestep for advection in v is dt/2 so we can choose dt_scalings larger than
% expected.
% Assumption on the transport vector, true for the considered tests: E is 
% always smaller than v (whose largest value is determined by the domain)
% 
% To understand the meaning of "CIRCULANT MATRIX" refer to the paper by
% Crouseilles et al. (2011)
%
% Polynomials degree for DG method: k; splines order: k+2 with highest
% regularity k+1 on internal and boundary knots (periodic boundary conditions)
%
% The problem we consider is: df/dt + v df/dx + E df/dv = f
% The source term is added only to the transport equation in x direction: it is
% given as a lambda function of x, v, dt, initial and final time for each
% timestep where we must compute it.
% 
% The code is structured as a nested loop: the outer one on the mesh
% refinements, the inner one on the dt scalings. The errors are computed at the
% end of the time loop for each scaling and saved in a vector; all the vectors
% are then collected in a struct associated to each refinement level.
%
% Change the fields fem, fem1D, nqn according to the test we want to perform.

data = struct('domain',      [-pi,pi; -4,4],...
              'test',        'manufactured',...
              'initial_f',   @(x,y) (2-cos(2*x)).*exp(-0.25*(4*y-1).^2),...
              'exact_sol',   @(x,y,t) (2-cos(2*x - 2*pi*t)).*exp(-0.25*(4*y-1).^2),...
              'source',      @(x,y,dt,tin,tf) 0.25*exp(-0.25*(4*y-1).^2)./(y-pi).*(((2*sqrt(pi)+1).*(4*y-2*sqrt(pi)).*(cos(2*x-2*y*dt-2*pi*tin)-cos(2*x-2*pi*tf))) + ...
                                                +(0.5*sqrt(pi)*(4*y-1)*(cos(2*x-2*pi*tf).^2-cos(2*x-2*y*dt-2*pi*tin).^2))),...
              'rho_0',       sqrt(pi),...   % Background density for Poisson
              'fem',         'Q2',...
              'fem1D',       'P2',...
              'BC',          'Period',...   % Boundary conditions in \Omega_x
              'nqn',         [3,3],...      % Number of quadrature nodes per direction (x,v)
              'time',         0,...
              'Tend',         1,...         % Final simulation time
              'uex_poisson',   @(x,t) cos(2*x-2*pi*t)*sqrt(pi)/8,...
              'gradex_poisson',@(x,t) -sqrt(pi)/4*sin(2*x-2*pi*t),...
              'type_mesh',    'CART');
          
nref = [3,4,5,6,7]; 
scalings = [2,4,8,16];  % It must be dt < dx/trasp 
[nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(data.nqn);
% Translate nodes on the reference element [0,1]. We can do it only once if we
% assume that we set the same number of nodes in both directions.
nodes_ref = (nodes_1D{1} + 1)*0.5; 

% Create input data structures to solve Poisson problem with IGA (GeoPDEs)
% Assumption: same 1D mesh and quandrature nodes used for DG method
xmin = data.domain(1,1); % Assuming the domain is a 1D interval 
xmax = data.domain(1,2);
pb_data.geo_name = nrbline ([xmin 0], [xmax 0]);
% Physical parameters
pb_data.c_diff  = @(x) ones(size(x));
% Average (of the solution) we want to impose (optional)
pb_data.avg     = 0; 

% REFINEMENTS LOOP
count = 0;
for ref = nref % Loop over refinements
    count = count + 1; % Needed for the struct to save errors
    fprintf(['\nSTART TEST ref = ' num2str(ref) '\n']);
              
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

    % MASS MATRIX 1D (inverse): local and global for the 2 directions
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
    velocity_values = flip(unique(femregionY.dof));
    temp = flip(femregionY.dof);
    for i = 1:length(velocity_values)
        unique_coord_dofY{i} = find( temp == velocity_values(i));
    end
    Vstruct = struct('values', velocity_values, 'unique_coord_dof', {unique_coord_dofY});
    
    % Helper structure to deal with E in building the matrix for the transport
    % in x direction. In the Vlasov-Poisson case we cannot precompute the 
    % values uX since they depend on time and on the previous solutions.
    dof1DX = unique(femregionX.dof);
    for i = 1:length(dof1DX)
        unique_coord_dofX{i} = find(femregionX.dof == dof1DX(i));
    end
    % To evaluate the IGA basis functions for E we must rescale the dofs
    % coordinates to the parametric domain (0,1)
    dof1DX_param = unique(dof_parametric_domain(femregionX));
    
    % ISOGEOMETRIC SPACE
    % MESH, SPACE, MATRIX FOR POISSON 1D (IGA code)
    method_data.degree     = femregion.degree +2;      % Degree of the splines
    method_data.regularity = femregion.degree +1;      % Regularity of the splines
    method_data.nsub       = 2^ref;                    % Number of subdivisions (same as DG methods)
    method_data.nquad      = data.nqn(1);              % Points for the Gaussian quadrature rule
    method_data.continuityBC = femregion.degree + 1;   % Continuity required at x=0 and x=1 (max by default)
    method_data.set_avg    = true;                     % Constraint on the solution to ensure uniqueness
    
    [geometry, msh_temp, space_temp, lap_mat] = matrix_laplace_periodicBC (pb_data, method_data);
    [L_lap, U_lap] = lu(lap_mat); % LU factorization
    msh = msh_precompute(msh_temp);
    space = sp_precompute(space_temp, msh, 'gradient',true);
    % basisIGAeval is a matrix with n_rows = n_basis_IGA and n_columns = n_dofDG_1D 
    [basisIGAeval, ~] = sp_eval_shape_functions(space_temp, geometry, {dof1DX_param}, {'gradient'});

    % CONSTRUCT SPACE AND MESH WITH TWICE THE NUMBER OF QUADRATURE NODES TO
    % COMPUTE ERRORS
    % Construct space with twice the number of quadrature nodes splines are
    % polynomial with higher order than DG, so even if they are integrated
    % exactly on the same number of qn, they may not be enough for the error
    % Perform degree elevation and knot insertion for the unclamping
    domain_elev = nrbdegelev(geometry.nurbs, method_data.degree  - 1); % number of times we need to increase the order (assuming we start from linear)
    [knots, zeta] = kntrefine (domain_elev.knots, 2^ref-1, method_data.degree, method_data.regularity);
    domain_clamped = nrbkntins(domain_elev, knots(method_data.degree+2 : end-method_data.degree-1));
    domain_unclamped = nrbunclamp(domain_clamped, method_data.continuityBC);

    % Construct msh and space structure with twice the number of quadrature nodes
    rule     = msh_gauss_nodes (2*data.nqn(1));
    [qn, qw] = msh_set_quad_nodes (zeta, rule);
    msh_err_temp      = msh_cartesian (zeta, qn, qw, geometry);
    space_err_temp    = sp_bspline (domain_unclamped.knots, method_data.degree, msh_err_temp); 
    msh_err = msh_precompute(msh_err_temp);
    space_err = sp_precompute(space_err_temp, msh_err, 'gradient',true);
        
    % Compute the initial ELECTRIC FIELD (phi = Poisson pb solution)
    % res is a matrix with n_rows = n_1Dnodes per element in x, and
    % n_columns = n_elements of the mesh in \Omega_x
    f_old = data.initial_f(femregion.dof(:,1), femregion.dof(:,2)); 
    [res, ~] = compute_vertical_line_integral(basis2D, f_old, femregion, region); 
    f_0 = data.rho_0*ones(size(res)) - res; % source term for Poisson
    phi_0 = solve_laplace1D_periodicBC_given_matrix(f_0, L_lap, U_lap, space, msh, method_data, pb_data);    
    vtk_pts = {linspace(0, 1, 100)}; 
    [E_0diagn, ~] = get_electric_field(space, msh, phi_0); % Here it is computed on the quadrature nodes, not on the dofs
    
%     [eu, F] = sp_eval (phi_0, space_temp, geometry, vtk_pts);
%     figure(); plot (F, eu,'r-'); title('Poisson solution at t = 0');
    
    % COMPUTE INITIAL AND CONSTANT VALUES FOR THE DIAGNOSTICS
    % We must also compute the mass matrix for the 2D DG space.
    M_dg2D = mass_matrix(femregion, basis2D);
    [energy0, charge0, L2norm0, ds, rhs_1, rhs_v2] = ...
            initialize_diagnostics(f_old, E_0diagn, M_dg2D, msh, femregion, basis2D );
    
    % Arrays storing errors for each dt scaling
    errL2vlasov = [];
    errL2poisson = [];
    errL2electric = [];
    times = [];
    dt_vector = [];
    
    k = 0;
    for dt_scaling = scalings % Loop over dt_scalings
        k = k+1;
        
        % Conserved properties: d_ = "relative difference with respect to initial value"
        d_energy = []; 
        d_charge = [];
        d_L2normVlasov = [];
        L2normE = []; % simple evolution of the L2 norm of the electric field
        
        % We take the smallest dx to find the dt satisfying our assumption on
        % the starting point for the characteristics.
        dxx = min(femregionX.h, femregionY.h); 
        dtt = dxx/dt_scaling;
        % Actually, we use the following smaller values (same as in the Eulerian method)
        % to perform sound comparisons)
        dx = region.h;
        dt = dx/dt_scaling;
        dt_vector = [dt_vector, dt];
        
        tic
        
        % BUILD GLOBAL MATRIX for the transport in X direction, collecting all
        % the systems for the 1D transport problems associated to the discrete
        % values of the transport coefficient.
        Av = global_matrixV(femregionX, femregionY, shape1D, nodes_ref, w_1D{1}, 0.5*dt, structures_circX, Vstruct); 
       
        % INITIALIZATION FOR TIME EVOLUTION
        f_old = data.initial_f(femregion.dof(:,1), femregion.dof(:,2));
        f_quarter = zeros(size(f_old));
        f_half = zeros(size(f_old));
        f_new = zeros(size(f_old));
    
        % Time loop
        for t = 0 : dt : data.Tend

            % COMPUTE SOURCES (half time step)
            s_half = @(x,y) data.source(x,y, 0.5*dt, t, t+0.5*dt);
            s_end = @(x,y) data.source(x,y, 0.5*dt, t+0.5*dt, t+dt);
            source_half = source1D(femregionX, s_half, basis1D, Vstruct);
            source_end = source1D(femregionX, s_end, basis1D, Vstruct);

            % ADVECTION IN X (HALF TIME STEP) 
            temp = MinvX * Av * f_old(connX) + MinvX * source_half; % f_old(connX) changes the order of the dofs, from 2D to the chosen 1D ordering 
            f_quarter(connX) = temp; % temp contains the values of the dofs ordered in 1D -> we want them in the 2D order

            % COMPUTE ELECTRIC FIELD from f_quarter 
            [res, ~] = compute_vertical_line_integral(basis2D, f_quarter, femregion, region);
            f_poisson = data.rho_0*ones(size(res)) - res; % source term for Poisson
            phi_temp = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);


%             [eu, F] = sp_eval (phi_temp, space_temp, geometry, vtk_pts);
%             figure(); plot (F, eu,'r-',F,dati.uex_poisson(F,t+0.5*dt),'--'); title(['Soluzione Poisson a t = ' num2str(t+0.5*dt)]);

            % BasisIGAeval is a matrix with n_rows = n_basis_IGA and n_columns = n_dofDG_1D 
            % storing the values of the basis functions on the coordinates
            % associated to the dofs of the DG space.
            E_half = basisIGAeval'*phi_temp;
            % BUILD structs storing the TRANSPORT VALUES
            Estruct = struct('values', E_half, 'unique_coord_dof', {unique_coord_dofX});

            % ADVECTION IN V (FULL TIME STEP)
            Ae = global_matrixE(femregionX, femregionY, shape1D, nodes_ref, w_1D{1}, dt, structures_circY, Estruct); 
            temp = MinvY * Ae * f_quarter(connY);
            f_half(connY) = temp;

            % ADVECTION IN X (HALF TIME STEP) 
            temp = MinvX * Av * f_half(connX) + MinvX * source_end;
            f_new(connX) = temp;

            f_old = f_new;

            % COMPUTE DIAGNOSTICS AT t+dt (relative values)
%             postprocessing(femregion, dati, f_new); 
%             pause()
            [E_half_diagn, ~] = get_electric_field(space, msh, phi_temp);
            [energy_cur, charge_cur, L2normVlas_cur, L2normE_cur] = compute_diagnostics(f_new, E_half_diagn, M_dg2D, ds, rhs_1, rhs_v2);
            d_energy = [d_energy, abs(energy_cur-energy0)/energy0];
            d_charge = [d_charge, (charge_cur-charge0)/charge0];
            d_L2normVlasov = [d_L2normVlasov, (L2normVlas_cur-L2norm0)/L2norm0];
            L2normE = [L2normE, L2normE_cur];
            
          
        end % Time loop

        diagnostics_scaling{k} = struct('type','base','d_energy_rel', d_energy, 'd_charge_rel', d_charge, ...
            'd_L2normVlasov_rel', d_L2normVlasov, 'L2normE', L2normE, 'scaling', dt_scaling, 'ref',ref);
         
        TT = toc;
        times = [times, TT];
        
        % COMPUTE ERRORS
        data.time = t+dt;
        err_vlasov = compute_errors(data, femregion, f_new);
        errL2vlasov = [errL2vlasov, err_vlasov.E_L2];

        [res,~] = compute_vertical_line_integral(basis2D, f_new, femregion, region); 
        f_poisson = data.rho_0*ones(size(res)) - res; % source term for Poisson
        phi_new = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);

        % We use twice the number of quadrature nodes to compute errors.
        [err_poisson, err_electric] = errorsL2_field_gradient(space_err_temp, msh_err_temp, phi_new, @(x) data.uex_poisson(x, t+dt), @(x) data.gradex_poisson(x, t+dt));
        errL2poisson = [errL2poisson, err_poisson];
        errL2electric = [errL2electric, err_electric];
        
        % POSTPROCESSING
        [eu, F] = sp_eval (phi_new, space_temp, geometry, vtk_pts);
        figure(); plot (F, eu,'r-', F, data.uex_poisson(F, t+dt),'b--'); 
        legend('Numeric','Exact')
        title(['Poisson solution at Tend, ref = ' num2str(ref)]);
               
    end % Loop over dt_scalings
    
    % Save errors
    diagnostics{count} = diagnostics_scaling;
    dt_ref{count} = struct('ref',ref, 'dt_values_scalings', dt_vector);
    errori_ref{count} = struct('ref',ref, 'h2D', femregion.h, 'h1D', femregionX.h, 'times', times ,...
        'vlasov_scaling', errL2vlasov, 'poisson_scaling', errL2poisson, 'electric_scaling', errL2electric);
    
end % Loop over refinements

fprintf('!!! END TEST !!!\n')

%% Convergence analysis
% (should be modified if more than one dt_scaling is tested)
figure()
subplot(1,3,1)
loglog( hh, errL2vlasov, 'o-', hh, hh.^(femregionX.degree - 1),'--')
legend('L2 norm',['h^' num2str(femregionX.degree - 1)]);
xlabel('h');
title(['Errors with polynomials ' femregionX.fem])

subplot(1,3,2)
loglog(hh, errL2poisson, '+-', hh, hh.^(femregionX.degree -1),'--')
xlabel('h');
legend('L2 norm Poisson',['h^' num2str(femregionX.degree -1)]);

subplot(1,3,3)
loglog(hh, errL2electric, '*-', hh, hh.^(femregionX.degree -1),'--')
xlabel('h');
legend('L2 norm electric',['h^' num2str(femregionX.degree -1)]);

pL2vlasov = log(errL2vlasov(2:end)./errL2vlasov(1:end-1))./log(hh(2:end)./hh(1:end-1))
pL2poisson = log(errL2poisson(2:end)./errL2poisson(1:end-1))./log(hh(2:end)./hh(1:end-1))
pL2electric = log(errL2electric(2:end)./errL2electric(1:end-1))./log(hh(2:end)./hh(1:end-1))

%% RESULTS ANALYSIS WITH DIFFERENT CFL CONDITIONS
% NB: we define the CFL number as: max(trasp)dt/dx = max(trasp)/dt_scaling
% where max(trasp) = 4 always since the electric field is always smaller in this
% test. We can consider CFL = 1 without breaking the assumption on the starting
% point of the characteristics, since the timestep we use for the advection in x
% (where the transport coefficient is larger) is dt/2.

% Ricostruisco matrice con errori per i vari scaling: nscaling x nref
% Construct matrix storing the errors associated to the scalings: 
% n_rows = n_scalings, n_columns = n_ref
vlasov = zeros(length(scalings), length(nref));
poisson = zeros(length(scalings), length(nref));
electric = zeros(length(scalings), length(nref));
hh2D = zeros(1,length(nref));
hh1D = zeros(1,length(nref));

for j = 1:length(nref)
    for  i = 1:length(scalings)
        vlasov(i,j) = errori_ref{j}.vlasov_scaling(i);
        poisson(i,j) = errori_ref{j}.poisson_scaling(i);
        electric(i,j) = errori_ref{j}.electric_scaling(i);
    end
    hh2D(j) = errori_ref{j}.h2D;
    hh1D(j) = errori_ref{j}.h1D;
end

% Compute relative errors
load('exactL2norms.mat')
vlasov = vlasov / L2normexactVlasov;
poisson = poisson / L2uex;
electric = electric / L2graduex;

%% PLOTS WITH THE ERRORS
colors = hsv(length(scalings));
markers = ['o','+','*','d']; % knowing we used 4 different dt_scalings
figure()

for i = 1:length(scalings)
    subplot(1,3,1)
    loglog(hh2D, vlasov(i,:),'-','Marker',markers(i),'Color',colors(i,:),'LineWidth',1,'MarkerSize',7)
    hold on; grid on
    leg1{i} = ['CFL = ', num2str(4/scalings(i),'%.2f')];
    
    subplot(1,3,2)
    loglog(hh2D, poisson(i,:),'-','Marker',markers(i),'Color',colors(i,:), 'LineWidth',1,'MarkerSize',7)
    hold on; grid on
    leg2{i} = ['CFL = ', num2str(4/scalings(i),'%.2f')];
    
    subplot(1,3,3)
    loglog(hh2D, electric(i,:),'-','Marker',markers(i),'Color',colors(i,:), 'LineWidth',1,'MarkerSize',7)
    hold on; grid on
    leg3{i} = ['CFL = ', num2str(4/scalings(i),'%.2f')];
end

if strcmp(data.fem, 'Q1')
    % Case Q1
    subplot(1,3,1); 
    loglog(hh2D, 100*hh2D.^(femregion.degree + 1),'k--','LineWidth',1); 
    leg1{i+1} = ['$O(h^' num2str(femregion.degree+1) ')$'];
    legend(leg1,'Interpreter','latex','Location','SE','FontSize',14); xlabel('h','FontSize',14); 
    ylabel('$\| f- f_h\|_{L^2} \, / \, \|f\|_{L^2}$','Interpreter','latex','FontSize',12);


    subplot(1,3,2); 
    loglog(hh2D, hh2D.^(femregion.degree+1),'k--','LineWidth',1); 
    leg2{i+1} = ['$O(h^' num2str(femregion.degree+1) ')$'];
    legend(leg2,'Interpreter','latex','Location','SE','FontSize',14); xlabel('h','FontSize',14); 
    ylabel('$\| \phi- \phi_h\|_{L^2} \, / \, \|\phi\|_{L^2}$','Interpreter','latex','FontSize',12)
    title(['\textbf{Polynomial degree r = ' num2str(femregion.degree) ' with splines degree ' num2str(femregion.degree+2) '}'],'Interpreter','latex','FontSize',18) 

    subplot(1,3,3);
    loglog(hh2D, hh2D.^(femregion.degree+1),'k--','LineWidth',1); 
    leg3{i+1} = ['$O(h^' num2str(femregion.degree+1) ')$'];
    legend(leg3,'Interpreter','latex', 'Location','SE','FontSize',14); xlabel('h','FontSize',14); 
    ylabel('$\| E- E_h\|_{L^2} \, / \, \|E\|_{L^2}$','Interpreter','latex','FontSize',12)

elseif strcmp(data.fem, 'Q2')
    % Case Q2
    if exist('dxx','var') % dt is the same as the Eulerian method with splitting
        subplot(1,3,1); 
        loglog(hh2D, 200*hh2D.^(femregion.degree + 1),'k--','LineWidth',1); 
        leg1{i+1} = ['$O(h^' num2str(femregion.degree+1) ')$'];
        legend(leg1,'Interpreter','latex','Location','SE','FontSize',14); xlabel('h','FontSize',14); 
        ylabel('$\| f- f_h\|_{L^2} \, / \, \|f\|_{L^2}$','Interpreter','latex','FontSize',12);

        subplot(1,3,2); 
        loglog(hh2D, 0.3*hh2D.^(femregion.degree),'k-.','LineWidth',1); 
        loglog(hh2D, 5*hh2D.^(femregion.degree + 2),'k--','LineWidth',1); 
        leg2{i+1} = ['$O(h^' num2str(femregion.degree) ')$'];
        leg2{i+2} = ['$O(h^' num2str(femregion.degree + 2) ')$'];
        legend(leg2,'Interpreter','latex','Location','SE','FontSize',14); xlabel('h','FontSize',14); 
        ylabel('$\| \phi- \phi_h\|_{L^2} \, / \, \|\phi\|_{L^2}$','Interpreter','latex','FontSize',12)
        title(['\textbf{Polynomial degree r = ' num2str(femregion.degree) ' with splines degree ' num2str(femregion.degree+2) '}'],'Interpreter','latex','FontSize',20) 

        subplot(1,3,3);
        loglog(hh2D, hh2D.^(femregion.degree),'k-.','LineWidth',1); 
        loglog(hh2D, 10*hh2D.^(femregion.degree+2),'k--','LineWidth',1); 
        leg3{i+1} = ['$O(h^' num2str(femregion.degree) ')$'];
        leg3{i+2} = ['$O(h^' num2str(femregion.degree+2) ')$'];
        legend(leg3,'Interpreter','latex', 'Location','SE','FontSize',14); xlabel('h','FontSize',14); 
        ylabel('$\| E- E_h\|_{L^2} \, / \, \|E\|_{L^2}$','Interpreter','latex','FontSize',12)

    else % dt is different of the one for the Eulerian method with splitting
        subplot(1,3,1); 
        loglog(hh2D, 100*hh2D.^(femregion.degree + 1),'k--','LineWidth',1); 
        leg1{i+1} = ['h^' num2str(femregion.degree+1)];
        legend(leg1,'Location','SE','FontSize',14); xlabel('h','FontSize',14); ylabel('Error L^2 VLASOV','FontSize',14);
        
        subplot(1,3,2); 
        loglog(hh2D, 0.01*hh2D.^(femregion.degree),'k--','LineWidth',1); 
        leg2{i+1} = ['h^' num2str(femregion.degree)];
        legend(leg2,'Location','SE','FontSize',14); xlabel('h','FontSize',14); ylabel('Error L^2 POISSON','FontSize',14);
        title(['Errors SEMILG-DG with polynomials ' data.fem],'FontSize',18) 
    
        subplot(1,3,3);
        loglog(hh2D, 0.01*hh2D.^(femregion.degree),'k--','LineWidth',1); 
        leg3{i+1} = ['h^' num2str(femregion.degree)];
        legend(leg3, 'Location','SE','FontSize',14); xlabel('h','FontSize',14); ylabel('Error L^2 ELECTRIC','FontSize',14);
    end
elseif strcmp(data.fem, 'Q3')
    if exist('dxx','var') % dt is the same as the Eulerian method with splitting
        subplot(1,3,1); 
        loglog(hh2D, 5*hh2D.^(femregion.degree),'k-.','LineWidth',1); 
        loglog( hh2D, 200*hh2D.^(femregion.degree+1 ),'k--','LineWidth',1);
        leg1{i+1} = ['$O(h^' num2str(femregion.degree) ')$'];
        leg1{i+2} = ['$O(h^' num2str(femregion.degree+1) ')$']; % se tolgo il commento sopra devo mettere i+2
        legend(leg1,'Interpreter','latex','Location','SE','FontSize',14); xlabel('h','FontSize',14); 
        ylabel('$\| f- f_h\|_{L^2} \, / \, \|f\|_{L^2}$','Interpreter','latex','FontSize',12);

        subplot(1,3,2); 
        loglog(hh2D, 0.2*hh2D.^(femregion.degree-1),'k-.','LineWidth',1); 
        leg2{i+1} = ['$O(h^' num2str(femregion.degree-1) ')$'];
        legend(leg2,'Interpreter','latex','Location','SE','FontSize',14); xlabel('h','FontSize',14); 
        ylabel('$\| \phi- \phi_h\|_{L^2} \, / \, \|\phi\|_{L^2}$','Interpreter','latex','FontSize',12)
        title(['\textbf{Polynomial degree r = ' num2str(femregion.degree) ' with splines degree ' num2str(femregion.degree+2) '}'],'Interpreter','latex','FontSize',20) 

        subplot(1,3,3);
        loglog(hh2D, 0.2*hh2D.^(femregion.degree-1),'k-.','LineWidth',1); 
        leg3{i+1} = ['$O(h^' num2str(femregion.degree-1) ')$'];
        legend(leg3, 'Interpreter','latex','Location','SE','FontSize',14); xlabel('h','FontSize',14); 
        ylabel('$\| E- E_h\|_{L^2} \, / \, \|E\|_{L^2}$','Interpreter','latex','FontSize',12)

    else % dt is different of the one for the Eulerian method with splitting
        subplot(1,3,1); 
        loglog(hh2D, 5*hh2D.^(femregion.degree-1),'k--', hh2D, 30*hh2D.^(femregion.degree ),'k-.','LineWidth',1); 
        leg1{i+1} = ['h^' num2str(femregion.degree-1)];
        leg1{i+2} = ['h^' num2str(femregion.degree)];
        legend(leg1,'Location','SE','FontSize',14); xlabel('h','FontSize',14); ylabel('Error L^2 VLASOV','FontSize',14);

        subplot(1,3,2); 
        loglog(hh2D, 0.05*hh2D.^(femregion.degree-1),'k--','LineWidth',1); 
        leg2{i+1} = ['h^' num2str(femregion.degree-1)];
        legend(leg2,'Location','SE','FontSize',14); xlabel('h','FontSize',14); ylabel('Error L^2 POISSON','FontSize',14);
        title(['Errors SEMILG-DG with polinomials ' data.fem],'FontSize',18) 
        
        subplot(1,3,3);
        loglog(hh2D, 0.05*hh2D.^(femregion.degree-1),'k--','LineWidth',1); 
        leg3{i+1} = ['h^' num2str(femregion.degree-1)];
        legend(leg3, 'Location','SE','FontSize',14); xlabel('h','FontSize',14); ylabel('Error L^2 ELECTRIC','FontSize',14);

    end
end

%% PLOTS WITH THE DIAGNOSTICS
colors = hsv(length(scalings));

% For a fixed refinement, plot all the st_scalings
for i = 1:length(nref)
    
    figure()
    for  j = 1:length(scalings)
               
        tt = linspace(0, data.Tend, length(diagnostics{i}{j}.d_energy_rel));
        subplot(2,2,1)
        semilogy(tt,abs(diagnostics{i}{j}.d_energy_rel),'Color', colors(j,:)); xlabel('t'); ylabel('log rel\_energy\_diff');
        hold on
        title(['Conserved properties: fem = ' data.fem ', spline deg = ' num2str(femregion.degree+2) ,...
            ', cells = ' num2str(2^diagnostics{i}{j}.ref) ],'FontSize',18);
        
        subplot(2,2,2)
        semilogy(tt,abs(diagnostics{i}{j}.d_charge_rel),'Color',colors(j,:)); xlabel('t'); ylabel('log rel\_charge\_diff');
        hold on

        subplot(2,2,3)
        semilogy(tt,abs(diagnostics{i}{j}.d_L2normVlasov_rel),'Color',colors(j,:)); xlabel('t'); ylabel('log rel\_L2norm\_diff');
        hold on

        subplot(2,2,4)
        semilogy(tt,abs(diagnostics{i}{j}.L2normE),'Color',colors(j,:)); xlabel('t'); ylabel('log L2norm\_E');
        hold on
        
        leg0{j} = ['scaling = ' num2str(diagnostics{i}{j}.scaling)];    
    end
    
    subplot(2,2,4);
    legend(leg0,'Location','best');
end

