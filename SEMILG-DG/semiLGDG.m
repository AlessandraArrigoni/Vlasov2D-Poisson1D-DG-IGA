% Function that implements the SEMILAGRANGIAN DG method with II order splitting
% (for 2D Vlasov equation) coupled with the IGA method (for 1D Poisson equation)
% or the LONG TIME SIMULATIONS, i.e. without exact solution. 
%
% Both the FEM and IGA spaces and matrices are precomputed at the beginnning 
% (if constant over time).
% 
% The problem we consider is: df/dt + v df/dx + E df/dv = 0
% We apply periodic in x and "compact support" (i.e. the solution
% is 0 outside the domain) in y boundary conditions.
%
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
% INPUT : struct (see input test files for the fields)
% OUTPUT : struct with one name per refinement level denoting the file where the
% results are saved; we save the time evolution of the diagnostics.

function [filenames] = semiLGDG(data)


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
    for ref = data.nref % Loop over refinements
        fprintf(['\nSTART TEST ref = ' num2str(ref) '\n']);
        
        clear ff phiphi 
        
        % Conserved properties: d_ = "relative difference with respect to initial value"
        d_energy = [];
        d_charge = [];
        d_L2normVlasov = [];
        L2normE = []; % simple evolution of the L2 norm of the electric field

        % Store output file name
        filename = [data.testname '_semiLGDG_' data.fem '_ref' num2str(ref) '.mat'];
        filenames.(['ref' num2str(ref)]) = filename;

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
        [rowsX, colsX, ~, valuesinvX] = mass_matrix1D_structure(femregionX, basis1D);
        MinvX = global_matrixMass(femregionX, rowsX, colsX, valuesinvX);

        [rowsY, colsY, ~, valuesinvY] = mass_matrix1D_structure(femregionY, basis1D);
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
            unique_coord_dofX{i} = find( femregionX.dof == dof1DX(i));
        end
        % To evaluate the IGA basis functions for E we must rescale the dofs
        % coordinates to the parametric domain (0,1)
        dof1DX_param = unique(dof_parametric_domain(femregionX));

        % ISOGEOMETRIC SPACE
        % MESH, SPACE, MATRIX FOR POISSON 1D (IGA code)
        method_data.degree     = data.splines;         % Degree of the splines
        method_data.regularity = data.splines - 1;     % Regularity of the splines
        method_data.nsub       = 2^ref;                % Number of subdivisions (same as DG method)
        method_data.nquad      = data.nqn(1);          % Points for the Gaussian quadrature rule
        method_data.continuityBC = data.splines - 1;   % Continuity required at x=0 and x=1 (max by default)
        method_data.set_avg    = true;   
       
        [geometry, msh_temp, space_temp, lap_mat] = matrix_laplace_periodicBC (pb_data, method_data);
        [L_lap, U_lap] = lu(lap_mat); % LU factorization
        msh = msh_precompute(msh_temp);
        space = sp_precompute(space_temp, msh, 'gradient',true);
        % basisIGAeval is a matrix with n_rows = n_basis_IGA and n_columns = n_dofDG_1D 
        [basisIGAeval, ~] = sp_eval_shape_functions(space_temp, geometry, {dof1DX_param}, {'gradient'});

        % Compute the initial ELECTRIC FIELD (phi = Poisson pb solution)
        % res is a matrix with n_rows = n_1Dnodes per element in x, and
        % n_columns = n_elements of the mesh in \Omega_x
        f_old = data.initial_f(femregion.dof(:,1), femregion.dof(:,2)); 
        f_quarter = zeros(size(f_old));
        f_half = zeros(size(f_old));
        f_new = zeros(size(f_old));
        [res, ~] = compute_vertical_line_integral(basis2D, f_old, femregion, region); 
        f_0 = data.rho_0*ones(size(res)) - res; % source term for Poisson
        phi_0 = solve_laplace1D_periodicBC_given_matrix(f_0, L_lap, U_lap, space, msh, method_data, pb_data);   
        [E_0diagn, ~] = get_electric_field(space, msh, phi_0); % Here it is computed on the quadrature nodes, not on the dofs
   
        % COMPUTE INITIAL AND CONSTANT VALUES FOR THE DIAGNOSTICS
        % We must also compute the mass matrix for the 2D DG space.
        M_dg2D = mass_matrix(femregion, basis2D);
        [energy0, charge0, L2norm0, ds, rhs_1, rhs_v2] = ...
                initialize_diagnostics(f_old, E_0diagn, M_dg2D, msh, femregion, basis2D );
            
        ff.first = f_old;
        phiphi.first = phi_0;
                                               
        % We take the smallest dx to find the dt satisfying our assumption on
        % the starting point for the characteristics.
        dxx = min(femregionX.h, femregionY.h); 
        dtt = dxx/data.dt_scaling;
        % Actually, we use the following smaller values (same as in the Eulerian method)
        % to perform sound comparisons)
        dx = region.h;
        dt = dx/data.dt_scaling;

        % BUILD GLOBAL MATRIX for the transport in X direction, collecting all
        % the systems for the 1D transport problems associated to the discrete
        % values of the transport coefficient.
        Av = global_matrixV(femregionX, femregionY, shape1D, nodes_ref, w_1D{1}, 0.5*dt, structures_circX, Vstruct); 
        
        % Save initial values to output file
        save(filename,'region','femregion','data','ref','dt',...
            'femregionX','femregionY','connX','connY','dxx','dx',...
            'method_data','pb_data','space','msh','geometry');
        
        % TIME LOOP
        count = 0;
        fprintf('--- Start time loop ---\n');
        
        for t = 0 : dt : data.Tend
            count = count + 1;
            
            tic
            
            % ADVECTION IN X (HALF TIME STEP) 
            temp = MinvX * Av * f_old(connX); % f_old(connX) changes the order of the dofs, from 2D to the chosen 1D ordering 
            f_quarter(connX) = temp; % temp contains the values of the dofs ordered in 1D -> we want them in the 2D order

            % COMPUTE ELECTRIC FIELD from f_quarter 
            [res, ~] = compute_vertical_line_integral(basis2D, f_quarter, femregion, region); 
            f_poisson = data.rho_0*ones(size(res)) - res; % source term for Poisson
            phi_temp = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);

            % BasisIGAeval is a matrix with n_rows = n_basis_IGA and n_columns = n_dofDG_1D 
            % storing the values of the basis functions on the coordinates
            % associated to the dofs of the DG space.
            E_half = basisIGAeval'*phi_temp;
            Estruct = struct('values', E_half, 'unique_coord_dof', {unique_coord_dofX});

            % ADVECTION IN V (FULL TIME STEP)
            Ae = global_matrixE(femregionX, femregionY, shape1D, nodes_ref, w_1D{1}, dt, structures_circY, Estruct); 
            temp = MinvY * Ae * f_quarter(connY);
            f_half(connY) = temp;

            % ADVECTION IN X (HALF TIME STEP) 
            temp = MinvX * Av * f_half(connX);
            f_new(connX) = temp;

            f_old = f_new;
            
            T(count) = toc;

            % COMPUTE DIAGNOSTICS AT t+dt (relative values)
            [E_half_diagn, ~] = get_electric_field(space, msh, phi_temp);
            [energy_cur, charge_cur, L2normVlas_cur, L2normE_cur] = compute_diagnostics(f_new, E_half_diagn, M_dg2D, ds, rhs_1, rhs_v2);
            d_energy = [d_energy, abs(energy_cur-energy0)/energy0];
            d_charge = [d_charge, (charge_cur-charge0)/charge0];
            d_L2normVlasov = [d_L2normVlasov, (L2normVlas_cur-L2norm0)/L2norm0];
            L2normE = [L2normE, L2normE_cur];

            % SAVE INTERMEDIATE SOLUTIONS to both problems
            % They are stored in a struct
            if(mod(count, data.damp_file/dt)==0)
                
                fprintf('\n-------- Start saving data for t = %d --------\n',t+dt);
                ff.(['n' num2str(count*dt/data.damp_file)]) = f_new; 
                phiphi.(['n' num2str(count*dt/data.damp_file)]) = phi_temp;
                save(filename, 'ff','phiphi','d_energy','d_charge','d_L2normVlasov','L2normE','-append');
                fprintf('-------- End saving data for t = %d --------\n',t+dt);
            end
            
        end % End time loop
 
        % Total time needed for the time loops (without time for saving and
        % compute the diagnostics)
        timeloopTime = sum(T); 
        
        % POSTPROCESSING
        diagnostics = struct('d_energy_rel', d_energy, 'd_charge_rel', d_charge, ...
            'd_L2normVlasov_rel', d_L2normVlasov,'L2normE', L2normE, 'ref', ref);

        ff.last = f_new;
        phiphi.last = phi_temp;
        save(filename, 'timeloopTime', 'diagnostics',...
            'ff','phiphi','-append');

    end % Loop refinements

    fprintf('\n--------END REFINEMENT LEVEL = %d --------\n',ref);

end