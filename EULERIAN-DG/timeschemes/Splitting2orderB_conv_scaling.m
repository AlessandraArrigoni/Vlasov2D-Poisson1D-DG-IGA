function [filenames] = Splitting2orderB_conv_scaling(Data)
% Function that implements the 2th order Crank Nicolson time scheme and the
% coupling between DG (for Vlasov equation) and IGA (for Poisson equation)
% methods for the CONVERGENCE TEST, i.e. when an exact solution is known and a
% source term for Vlasov equation must be computed.
% Different values of dt_scaling (see input data file) are tested for each
% refinement.
%
% A second order SPLITTING technique is employed, where the first problem is the
% advection in v-direction, i.e. the equation:
% df / dt + E * df / dv = f
%
% Both the FEM and IGA spaces and matrices are precomputed at the beginnning 
% (if constant over time).
%
% INPUT: struct (see input test files for the fields Must contain a field "scalings")
% OUTPUT: struct with one name per refinement level denoting the file where the
% results are saved; we save the time evolution of the diagnostics and the error
% values at the final time. 

    % Create input data structures to solve Poisson problem with IGA (GeoPDEs)
    % Assumption: same 1D mesh and quandrature nodes used for DG method
    xmin = Data.domain(1,1); % Assuming the domain is a 1D interval 
    xmax = Data.domain(1,2);
    pb_data.geo_name = nrbline ([xmin 0], [xmax 0]);
    % Physical parameters
    pb_data.c_diff  = @(x) ones(size(x));
    % Average (of the solution) we want to impose (optional)
    pb_data.avg     = 0; 
    
    for ref = Data.nref % Loop over refinements
        clear ff phiphi % otherwise the values are added to the structs of the previous refinement

        % Conserved properties: d_ = "relative difference wrt initial value"
        d_energy = [];
        d_charge = [];
        d_L2normVlasov = [];
        L2normE = []; % simple evolution of the L2 norm of the electric field

        % Store output file name
        filename = [Data.testname '_splittingB_' Data.fem '_ref' num2str(ref) '.mat'];
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
        M_dg = mass_matrix(femregion, basis);
        
        
        % MESH, SPACE, MATRIX FOR POISSON 1D (IGA code)
        method_data.degree     = Data.splines;         % Degree of the splines
        method_data.regularity = Data.splines - 1;     % Regularity of the splines
        method_data.nsub       = 2^ref;                % Number of subdivisions (same as DG method)
        method_data.nquad      = Data.nqn(1);          % Points for the Gaussian quadrature rule
        method_data.continuityBC = Data.splines - 1;   % Continuity required at x=0 and x=1 (max by default)
        method_data.set_avg    = true;                 % Constraint on the solution to ensure uniqueness

        [geometry, msh_temp, space_temp, lap_mat] = matrix_laplace_periodicBC (pb_data, method_data);
        [L_lap, U_lap] = lu(lap_mat); % LU Factorization
        msh = msh_precompute(msh_temp);
        space = sp_precompute(space_temp, msh, 'gradient',true);
        
        % Save initial values to output file
        save(filename,'region','femregion','Data','ref',...
            'method_data','pb_data','space','msh','geometry');

        % Vectors to store errors at final time associated to the different
        % dt_scalings we test.
        errVlasov = []; 
        errPoisson = [];
        errElectric = [];
        
        k = 0;
        for dt_scaling = Data.scalings % Loop over dt_scalings
            k = k+1;
            fprintf('dt scaling = %d --------\n',dt_scaling);
            
            dt = region.h/dt_scaling; % time discretization parameter

            % CONSTANT PART OF THE TRANSPORT MATRIX
            Lv = matrix_onlyV(femregion,neighbour,Data,basis);        
            LHS2 = M_dg + 0.5*dt*Lv; 
            % LU factorization of the operator for the second step in the
            % splitting technique
            [L2, U2] = lu(LHS2); 
            
            % INITIALIZATION FOR TIME EVOLUTION
            f0 = Data.exact_sol(femregion.dof(:,1),femregion.dof(:,2),0);
            f_old = f0;
            f_new = zeros(femregion.ndof, 1);
            f_quarter = zeros(femregion.ndof, 1);
            f_half = zeros(femregion.ndof, 1);
            % Compute the initial ELECTRIC FIELD (phi = Poisson pb solution)
            % res is a matrix with n_rows = n_1Dnodes per element in x, and
            % n_columns = n_elements of the mesh in \Omega_x
            [res, ~] = compute_vertical_line_integral(basis, f0, femregion, region);
            f_poisson = Data.rho_0*ones(size(res)) - res; % source term for Poisson
            phi0 = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);
            [E_old, ~] = get_electric_field(space, msh, phi0); % Initial electric field

            % Compute the term depending on E for Vlasov equation (only to start the
            % loop, then we will use the matrix from step 3 in the algorithm)
            Le = matrix_onlyE(femregion,neighbour,Data,basis,E_old) ;
            LHS1 = M_dg + 0.25*dt*Le;
    
            ff.first = f0;
            phiphi.first = phi0;

            % COMPUTE INITIAL AND CONSTANT VALUES FOR THE DIAGNOSTICS
            [energy0, charge0, L2norm0, ds, rhs_1, rhs_v2] = ...
                initialize_diagnostics(f0, E_old, M_dg, msh, femregion, basis );
        
        
            count  = 0; T = [];
            for t = 0 : dt : Data.Tend-dt % Time loop
                count = count + 1;

                tic

                source0 = source_rhs_given_f(femregion, @(x,y) Data.source(x,y,t),basis);
                source1 = source_rhs_given_f(femregion, @(x,y) Data.source(x,y,t+dt),basis);

                % STEP 1: advection in v (beta = [0, E]), WITHOUT source term   
                rhs1 = (M_dg - 0.25*dt*Le)*f_old ;

                f_quarter = LHS1\rhs1;

                % STEP 2: advection in x (beta = [v, 0]), WITH source term
                rhs2 = (M_dg - 0.5*dt*Lv)*f_quarter + 0.5*dt*(source0 + source1);

                ytemp = L2\rhs2;
                f_half = U2\ytemp;

                % Compute E(n+1) from f_half
                [res, ~] = compute_vertical_line_integral(basis, f_half, femregion, region); 
                f_poisson = Data.rho_0*ones(size(res)) - res; % source term for Poisson
                phi_new = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);
                [E_new, ~] = get_electric_field(space, msh, phi_new); 

                % STEP 3: advection in v (beta = [0, E]), WITHOUT source term  
                % Compute the term depending on E for Vlasov equation
                Le = matrix_onlyE(femregion,neighbour,Data,basis,E_new) ;

                LHS3 = M_dg + 0.25*dt*Le;
                rhs3 = (M_dg - 0.25*dt*Le)*f_half;

                f_new = LHS3\rhs3;

                % Update solution.
                f_old = f_new;
                LHS1 = LHS3;

                % Save one iteration's time
                T(count) = toc;

                % COMPUTE DIAGNOSTICS (relative values)
                [energy_cur, charge_cur, L2normVlas_cur] = compute_diagnostics(f_new, E_new, M_dg, ds, rhs_1, rhs_v2);
                d_energy = [d_energy, abs(energy_cur-energy0)/energy0];
                d_charge = [d_charge, (charge_cur-charge0)/charge0];
                d_L2normVlasov = [d_L2normVlasov, (L2normVlas_cur-L2norm0)/L2norm0];

                % Save intermediate solutions to both problems
                % They are stored in a struct
                if(mod(count, 0.25/dt)==0)
                    % Compute solution of Poisson problem at the end of the timestep
                    [res, ~] = compute_vertical_line_integral(basis, f_new, femregion, region); 
                    f_poisson = Data.rho_0*ones(size(res)) - res; % source term for Poisson
                    phi_new = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);

                    ff.(['n' num2str(count)]) = f_new; 
                    phiphi.(['n' num2str(count)]) = phi_new;
                    save(filename, 'ff','phiphi','d_energy','d_charge','d_L2normVlasov','-append');
                end

            end % End time loop
            
            % Total time needed for the time loops (without time for saving and
            % compute the diagnostics)
            timeloopTime{k} = sum(T); 
            
            % POSTPROCESSING: output is a cell array with n_dt_scalings
            % elements; each element is a struct with the diagnostics.
            diagnostics{k} = struct('d_energy_rel', d_energy, 'd_charge_rel', d_charge, ...
                'd_L2normVlasov_rel', d_L2normVlasov, 'ref', ref, 'scaling',dt_scaling);

            % Find electric potential and field at the end of the simulation to
            % compute the errors in the most accurate way.
            [res, ~] = compute_vertical_line_integral(basis, f_new, femregion, region); 
            f_poisson = Data.rho_0*ones(size(res)) - res; % source term for Poisson
            phi_new = solve_laplace1D_periodicBC_given_matrix(f_poisson, L_lap, U_lap, space, msh, method_data, pb_data);

            ff.last = f_new;
            phiphi.last = phi_new;
            save(filename, 'timeloopTime', 'diagnostics',...
                'ff','phiphi','-append');
        
            % ERRORS: output is a vector for each type of error, storing the
            % errors at final time for all the considered dt_scaling for the
            % current refinement level.
            if Data.computeErrors
                Data.time = t+dt;
                err_vlasov = compute_errors(Data, femregion, f_new);
                errVlasov = [errVlasov, err_vlasov.E_L2];

                [poisson_l2, electric_l2] = errorsL2_field_gradient (space_temp, msh_temp, phi_new, @(x) Data.uex_poisson(x,t+dt), @(x) Data.gradex_poisson(x,t+dt));

                errPoisson = [errPoisson, poisson_l2];
                errElectric = [errElectric, electric_l2];
                
            end
            
        end % End loop over dt_scalings
        
        save(filename, 'errVlasov','errPoisson','errElectric','-append'); 
        fprintf('\n--------END REFINEMENT LEVEL = %d --------\n',ref);
        
    end % End loop over refinements


end

