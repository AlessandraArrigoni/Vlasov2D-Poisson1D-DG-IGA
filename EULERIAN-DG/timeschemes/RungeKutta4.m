function [filenames] = RungeKutta4(Data)
% Function that implements the 4th order Runge-Kutta time scheme and the
% coupling between DG (for Vlasov equation) and IGA (for Poisson equation)
% methods for the LONG TIME SIMULATIONS, i.e. without exact solution. 
%
% Both the FEM and IGA spaces and matrices are precomputed at the beginnning 
% (if constant over time).
%
% INPUT: struct (see input test files for the fields)
% OUTPUT: struct with one name per refinement level denoting the file where the
% results are saved; we save the time evolution of the diagnostics.

    % Create input data structures to solve Poisson problem with IGA (GeoPDEs)
    % Assumption: same 1D mesh and quandrature nodes used for DG method
    xmin = Data.domain(1,1); % Assuming the domain is a 1D interval 
    xmax = Data.domain(1,2);
    pb_data.geo_name = nrbline ([xmin 0], [xmax 0]);
    % Physical parameters
    pb_data.c_diff  = @(x) ones(size(x));
    % Average (of the solution) we want to impose (optional)
    pb_data.avg     = 0; 
    
    for ref = Data.nref

        % Conserved properties: d_ = "relative difference with respect to initial value"
        d_energy = [];
        d_charge = [];
        d_L2normVlasov = [];
        L2normE = []; % simple evolution of the L2 norm of the electric field

        % Store output file name
        filename = [Data.testname '_rk4_' Data.fem '_ref' num2str(ref) '.mat'];
        filenames.(['ref' num2str(ref)]) = filename;
        
        fprintf('\n--------REFINEMENT LEVEL = %d --------\n',ref);

        % MESH AND SPACE FOR VLASOV 2D (DG method)
        region = generate_mesh(Data, ref);
        femregion = create_dof(Data, region);
        if ref < 4
            neighbour = neighbours_Period(femregion); 
        else % load saved neighbour structure to save computational time
            load(['neighbour_ref' num2str(ref) '.mat']);
        end

        [shape_basis] = basis_lagrange(femregion.fem); % scalar shape functions
        [nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(Data.nqn);
        [dphiq, Grad, B_edge] = evalshape(femregion,shape_basis,nodes_2D,nodes_1D,femregion.nln);


        % PERMUTE MATRICES TO ALLOW FASTER COMPUTATIONS
        gradphi = permute(Grad,[2,3,1]); % 2 x nln x nqn2D
        dphiq = reshape(dphiq,[length(nodes_2D),femregion.nln]); % nqn2D x nln

        basis = struct('nodes_1D', {nodes_1D}, 'w_1D', {w_1D}, 'nodes_2D', nodes_2D, 'w_2D', w_2D, ...
                'dphiq', dphiq, 'gradphi', gradphi, 'B_edge', {B_edge});

        % MASS MATRIX FOR THE TIME EVOLUTION
        [M_dg, Minv] = mass_matrix(femregion, basis);

        % CONSTANT PART OF THE TRANSPORT MATRIX
        Lv = matrix_onlyV(femregion, neighbour, Data, basis); 

        % MESH, SPACE, MATRIX FOR POISSON 1D (IGA code)
        method_data.degree     = Data.splines;         % Degree of the splines
        method_data.regularity = Data.splines - 1;     % Regularity of the splines
        method_data.nsub       = 2^ref;                % Number of subdivisions (same of DG method)
        method_data.nquad      = Data.nqn(1);          % Points for the Gaussian quadrature rule
        method_data.continuityBC = Data.splines - 1;   % Continuity required at x=0 and x=1 (max by default)
        method_data.set_avg    = true;                 % Constraint on the solution to ensure uniqueness

        [geometry, msh_temp, space_temp, lap_mat] = matrix_laplace_periodicBC (pb_data, method_data);
        [L_lap, U_lap] = lu(lap_mat); % LU factorization
        msh = msh_precompute(msh_temp);
        space = sp_precompute(space_temp, msh, 'gradient',true);

        % INITIALIZATION FOR TIME EVOLUTION
        f0 = Data.initial_f(femregion.dof(:,1),femregion.dof(:,2));
        f_old = f0;
        f_new = zeros(femregion.ndof, 1);
        kappa1 = zeros(femregion.ndof, 1);
        kappa2 = zeros(femregion.ndof, 1);
        kappa3 = zeros(femregion.ndof, 1);
        kappa4 = zeros(femregion.ndof, 1);
        % Compute the initial ELECTRIC FIELD (phi = Poisson pb solution)
        % res is a matrix with n_rows = n_1Dnodes per element in x, and
        % n_columns = n_elements of the mesh in \Omega_x
        [res, ~] = compute_vertical_line_integral(basis, f0, femregion, region); 
        f_poisson = Data.rho_0*ones(size(res)) - res; % source term for Poisson
        phi0 = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);
        [E_old, ~] = get_electric_field(space, msh, phi0); % Initial electric field

        % COMPUTE INITIAL AND CONSTANT VALUES FOR THE DIAGNOSTICS
        [energy0, charge0, L2norm0, ds, rhs_1, rhs_v2] = ...
            initialize_diagnostics(f0, E_old, M_dg, msh, femregion, basis );

        dt = region.h/Data.dt_scaling; % time discretization parameter (h/4 per Q2, h/6 per Q3 for stability reasons)

        % Save initial values to output file
        save(filename,'region','femregion','Data','ref',...
            'method_data','pb_data','space','msh','geometry','f0','phi0','dt');

        count  = 0;

        for t = 0 : dt : Data.Tend-dt % Time loop
            count = count + 1;

            tic

            fprintf('\n--------STEP 1 for t = %d \n',t);
            % Compute the term depending on E for Vlasov equation
            Le = matrix_onlyE(femregion, neighbour, Data ,basis, E_old);       
            LDG1 = Lv + Le;

            % FIRST step in Runge Kutta 4
            kappa1 = dt*Minv*(-LDG1*f_old); 

            % Compute source for Poisson at each intermidiate step from the
            % solution of Vlasov equation:
            % compute E(n+1/4) from da 0.5*kappa1 + f_old
            [res, ~] = compute_vertical_line_integral(basis, 0.5*kappa1 + f_old, femregion, region); 
            f_poisson = Data.rho_0*ones(size(res)) - res; % source term for Poisson
            phi_temp = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);
            [E_temp, ~] = get_electric_field(space, msh, phi_temp);


            fprintf('--------STEP 2 for t = %d \n',t);
            % Compute the term depending on E for Vlasov equation
            Le = matrix_onlyE(femregion, neighbour, Data ,basis, E_temp);   
            LDG2 = Lv + Le;

            % SECOND step in Runge Kutta 4
            kappa2 = dt*Minv*(-LDG2*(f_old + 0.5*kappa1));

            % Compute E(n+1/2) from 0.5*kappa2 + f_old
            [res, ~] = compute_vertical_line_integral(basis, 0.5*kappa2 + f_old, femregion, region);  
            f_poisson = Data.rho_0*ones(size(res)) - res; % source term for Poisson
            phi_temp = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);
            [E_temp, ~] = get_electric_field(space, msh, phi_temp); 

            fprintf('--------STEP 3 for t = %d \n',t);
            % Compute the term depending on E for Vlasov equation
            Le = matrix_onlyE(femregion, neighbour, Data ,basis, E_temp);   
            LDG3 = Le + Lv;

            % THIRD step in Runge Kutta 4
            kappa3 = dt*Minv*(-LDG3*(f_old + 0.5*kappa2));

            % Compute E(n+3/4) from kappa3 + f_old
            [res, ~] = compute_vertical_line_integral(basis, kappa3 + f_old, femregion, region); 
            f_poisson = Data.rho_0*ones(size(res)) - res; % source term for Poisson
            phi_temp = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);
            [E_temp, ~] = get_electric_field(space, msh, phi_temp);

            fprintf('--------STEP 4 for t = %d \n',t);
            % Compute the term depending on E for Vlasov equation
            Le = matrix_onlyE(femregion, neighbour, Data ,basis, E_temp);   
            LDG4 = Lv + Le;

            % FOURTH step in Runge Kutta 4
            kappa4 = dt*Minv*(-LDG4*(f_old + kappa3));

            % Update solution: linear combination of intermediate steps
            f_new = f_old + 1/6*(kappa1 + 2*kappa2 + 2*kappa3 + kappa4);
            f_old = f_new;

            % Compute electric field at step n+1
            [res, ~] = compute_vertical_line_integral(basis, f_new, femregion, region);  
            f_poisson = Data.rho_0*ones(size(res)) - res; % source term for Poisson
            phi_new = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);
            [E_new, ~] = get_electric_field(space, msh, phi_new); 
            % Transport field for Vlasov equation
            E_old = E_new;

            % Save one iteration's time
            T(count) = toc;

            % COMPUTE DIAGNOSTICS (relative values)
            [energy_cur, charge_cur, L2normVlas_cur, L2normE_cur] = compute_diagnostics(f_new, E_new, M_dg, ds, rhs_1, rhs_v2);
            d_energy = [d_energy, abs(energy_cur-energy0)/energy0];
            d_charge = [d_charge, (charge_cur-charge0)/charge0];
            d_L2normVlasov = [d_L2normVlasov, (L2normVlas_cur-L2norm0)/L2norm0];
            L2normE = [L2normE, L2normE_cur];

            % SAVE INTERMEDIATE SOLUTIONS to both problems
            % They are stored in a struct
            if(mod(count, Data.damp_file/dt)==0)
                ff.(['n' num2str(count)]) = f_new; 
                EE.(['n' num2str(count)]) = E_new;
                phiphi.(['n' num2str(count)]) = phi_new;
                save(filename, 'ff','EE','phiphi','L2normE','d_energy','d_charge','d_L2normVlasov','-append');
            end

        end % End time loop

        % Total time needed for the time loops (without time for saving and
        % compute the diagnostics)
        timeloopTime = sum(T); 
        
        % POSTPROCESSING
        diagnostics = struct('d_energy_rel', d_energy, 'd_charge_rel', d_charge, ...
            'd_L2normVlasov_rel', d_L2normVlasov, 'L2normE', L2normE, 'ref', ref);

        ff.last = f_new;
        EE.last = E_new;
        phiphi.last = phi_new;
        save(filename, 'timeloopTime', 'diagnostics',...
            'ff','EE','phiphi','-append');
        
        fprintf('\n--------END REFINEMENT LEVEL = %d --------\n',ref);
    end

end