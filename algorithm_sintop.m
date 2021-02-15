% BENCHMARK CASE: always a sin() on the top boundary.

insert(py.sys.path,int32(0),'/Users/alearrigoni/defeaturing'); % update python path
% % HOW TO RELOAD THE MODULE inout.py IF I CHANGE SOMETHING (not sure it
% works)
% clear classes;
% mod = py.importlib.import_module('inout');
% py.importlib.reload(mod);


% INITIALIZE PROBLEM and cost functional and choose parameters
max_iter_algo = 15;         % Defeaturing algorithm 
defeaturing_tol = 1e-6;     % Tolerance on the shape derivatives estimator
exactgeom_thresh = 1e-8;    % Tolerance on the boundary approximation for the exact geometry
proj_type = 'H1';           % or 'L2'. How to approximate the boundary
sin_parameters = [0.1, 7, 3]; % amplitude, frequency, width (see Jochen's code)
Jstring = 'u_squared';      % 'u','u_squared','grad_usquared' or...
Jdomain = '0';              % 0 = interior, 1,2,3,4 = left, right, bottom, top boundary
est_shape_der = [];
Jvalues = [];
ndofs_geom = [];

% PHYSICAL DATA OF THE PROBLEM: source = 0 and we set an inflow boundary on the left 
clear problem_data  

% Type of boundary conditions for each side of the domain: we need to pass
% them to python for the shape derivative!
problem_data.nmnn_sides   = [1];
problem_data.drchlt_sides = [2 3 4];
boundary_cond = {['N', '1'], ['D','0'], ['D','0'], ['D','0']};

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));
problem_data.grad_c_diff = @(x, y) cat (1, ...
            reshape (zeros(size(x)), [1, size(x)]), ...
            reshape (zeros(size(x)), [1, size(x)]));
        
% Source and boundary terms
problem_data.f = @(x, y) zeros (size (x));
problem_data.g = @neumann_BC; % Inflow boundary on the left
problem_data.h = @(x, y, ind) 0.*x.*y;
% Exact solution (optional)
problem_data.uex     = [];

% ADJOINT PROBLEM DATA (depend on the data for the pde and on the functional)
% NB: pay attention to the sign of f and g, since we will solve -div()=f!
clear adjproblem_data
adjproblem_data = get_adjoint_problem_data(problem_data, Jstring, Jdomain);

% METHOD PARAMETERS
clear method_data
method_data.degree      = [3 3];       % Degree of the splines
method_data.regularity  = [2 2];       % Regularity of the splines
method_data.nsub_coarse = [2 2];       % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];       % Number of subdivisions for each refinement: always dyadic
method_data.nquad       = [4 4];       % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';  % 'simplified' (only children functions) or 'standard' (full basis): always 'standard'
method_data.truncated   = 1;           % 0: False, 1: True. Always 1 here.

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
% adaptivity_data.flag = 'functions';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = .5;
adaptivity_data.mark_strategy = 'GERS'; % GR = uniform refinement
adaptivity_data.max_level = 10;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 5;
adaptivity_data.max_nel = 15000;
adaptivity_data.tol = 1e-6;

% CREATE COARSEST MESH of the parametric domain and the corresponding
% space, to start the algorithm.
nel_lev0 = method_data.nsub_coarse;
sp_degree = method_data.degree;
active_elements_geom = {[1:prod(method_data.nsub_coarse)]'};

% Only needed to use other functions of geopdes
geom_base = geo_load ('geo_square.txt');
[knots, zeta] = kntrefine (geom_base.nurbs.knots, method_data.nsub_coarse-1, method_data.degree, method_data.regularity);
  
rule     = msh_gauss_nodes (method_data.nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh_base   = msh_cartesian (zeta, qn, qw, geom_base);
space_base = sp_bspline (knots, method_data.degree, msh_base);
hmsh_base  = hierarchical_mesh (msh_base, method_data.nsub_refine); 
hspace_base = hierarchical_space (hmsh_base, space_base, method_data.space_type, method_data.truncated);

% Compute ONCE (and then store and load!) the "exact" geometry and the
% "exact" numerical solution for the "exact" value of the functional
filename_exact = 'exact_geometries/exact_geom8.mat';
[Jexact, filename_exact_out] = exact_geometry_and_functional(filename_exact, exactgeom_thresh, nel_lev0, sp_degree, proj_type, sin_parameters, ...
                                                problem_data, adaptivity_data, Jstring, Jdomain, hmsh_base, hspace_base, false);

% DEFEATURING LOOP
for iter = 1:max_iter_algo
    
    filename_geom = ['inoutfiles/geometry_iter' num2str(iter) '.mat'];
    save(filename_geom, 'active_elements_geom','nel_lev0','sp_degree','proj_type','sin_parameters');
    disp([newline 'COMPUTING PARAMETRIZATION // ITER ' num2str(iter)]);
    py.inout.approximate_geometry(filename_geom) % call python function
    load(['inoutfiles/geometry_iter' num2str(iter) '_out.mat']);
    if ~iscell(active_elements_geom_nofold)
        active_elements_geom_nofold = {active_elements_geom_nofold};
    end
       
    % Create initial mesh and space (corresponding to the geometry).
    % Always start from the base case since the geometry approximation
    % might also require the refinement of intermediate levels, not only of
    % the last one from the previous iteration.
    nlevels = numel(active_elements_geom_nofold)-1; % active cells in python output file
    adaptivity_data_geom.flag = 'elements';
    hmsh_geom = hmsh_base; hspace_geom = hspace_base;
    for ref = 1:nlevels
        marked_geom = cell (ref,1);
        marked_geom{ref} = setdiff(hmsh_geom.active{hmsh_geom.nlevels}, active_elements_geom_nofold{ref}); % it also sorts the vector 
        [hmsh_geom, hspace_geom] = adaptivity_refine (hmsh_geom, hspace_geom, marked_geom, adaptivity_data_geom); 
    end

    % GEOMETRY RECONSTRUCTION: define a new geometry with nurbs toolbox.
    % Write the control points associated to the basis of the finest level
    % (tensor product space on the unit square) and contruct a geometry.
    % Save the mapping in the initial hmesh used to solve the pde.
    C = hspace_subdivision_matrix (hspace_geom);
    fine_cpoints =  C{hspace_geom.nlevels} * cpoints; 
    fine_cpoints = reshape(fine_cpoints', [hmsh_geom.rdim, hspace_geom.space_of_level(end).ndof_dir]);
    % Build the geometry mappings
    fine_knots = hspace_geom.space_of_level(end).knots;
    nurbs = nrbmak(fine_cpoints, fine_knots); 
    geometry_iter = geo_load (nurbs);    
    % Update the mapping
    hmsh = hmsh_geom; hspace = hspace_geom;
    hmsh = update_mappings(hmsh, geometry_iter);
    
    % Visualize mesh and geometry
%     hmsh_plot_cells(hmsh_geom,10,figure); view(2);
%     title(['Geometric mesh no fold ITER ' num2str(iter)],'FontSize',16);
    hmsh_plot_cells(hmsh,10,figure); view(2);
    title(['Geometry ITER ' num2str(iter)],'FontSize',16);
    
    
        % ADAPTIVE LOOP: compute u and p (dual solution) on the same mesh.
        it_adapt = 0; gest = [];
        while true
            it_adapt = it_adapt + 1;    
            if (~hspace_check_partition_of_unity (hspace, hmsh))
                disp('ERROR: The partition-of-the-unity property does not hold.')
                break
            end
            
            % Solve PDE
            [sol_u, lap_mat] = adaptivity_solve_laplace (hmsh, hspace, problem_data);
            nel(it_adapt) = hmsh.nel; ndof(it_adapt) = hspace.ndof;
            % Estimate and check
            est = adaptivity_estimate_laplace (sol_u, hmsh, hspace, problem_data, adaptivity_data);
            gest(it_adapt) = norm (est); % save values of the estimator
            if (gest(it_adapt) < adaptivity_data.tol)
                disp('Success: The error estimation reached the desired tolerance'); 
                break
            elseif (it_adapt == adaptivity_data.num_max_iter)
                disp('Warning: reached the maximum number of iterations')
                break
            elseif (hmsh.nlevels >= adaptivity_data.max_level)
                disp('Warning: reached the maximum number of levels')
                break
            elseif (hspace.ndof > adaptivity_data.max_ndof)
                disp('Warning: reached the maximum number of DOFs')
                break
            elseif (hmsh.nel > adaptivity_data.max_nel)
                disp('Warning: reached the maximum number of elements')
                break
            end
            % Mark and refine
            [marked, ~] = adaptivity_mark (est, hmsh, hspace, adaptivity_data);            
            [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
        end
        
    % Save overkilled mesh and space (clearer, for now)
    hmsh_ref = hmsh; 
    hspace_ref = hspace;

    % Solve adjoint problem (on the same mesh and space as the pde, for now)
    sol_p = solve_adjoint_laplace(hmsh_ref, hspace_ref, lap_mat, sol_u, adjproblem_data);
        
    % COMPUTE FUNCTIONAL
    J = compute_functionalJ(Jstring, Jdomain, sol_u, hmsh_ref, hspace_ref);
    Jvalues = [Jvalues, J];
    ndofs_geom = [ndofs_geom, hspace_geom.ndof];
    
%     % Visualize solution and adjoint solution
%       figure(); subplot(1,2,1); title(['Solutions u and p ITER ' num2str(iter)],'FontSize',16);
%       sp_plot_solution(sol_u, hspace_ref, geometry_iter, 101); view(2); colorbar();
%       subplot(1,2,2);
%       sp_plot_solution(sol_p, hspace_ref, geometry_iter, 101); view(2); colorbar();

%     % Save and export solutions to vtk
%     filename_u = ['plots/sol_iter' num2str(iter)];
%     filename_p = ['plots/sol_dual_iter' num2str(iter)];
%     sp_to_vtk(sol_u, hspace_ref, geometry_iter, 151, filename_u, 'u');
%     sp_to_vtk(sol_p, hspace_ref, geometry_iter, 151, filename_p, 'p');
       
    % COMPUTE SHAPE DERIVATIVES: vector of hspace_geom_ref.ndof x 2
    active_elements_ref = hmsh_ref.active; 
    filename_der = ['inoutfiles/shapeder_iter' num2str(iter) '.mat'];
    save(filename_der,'sol_u','sol_p','Jstring','Jdomain','boundary_cond',...
        'active_elements_geom','active_elements_geom_nofold','active_elements_ref','cpoints',...
        'nel_lev0', 'sp_degree','proj_type','sin_parameters');
    disp([newline 'COMPUTING SHAPE DERIVATIVES // ITER ' num2str(iter)]);
    py.inout.shape_derivative(filename_der); % call python function
    load(['inoutfiles/shapeder_iter' num2str(iter) '_out.mat']);
    
    % CHECK ESTIMATOR: if smaller than tolerance, stop.
    % For now we consider the MAXIMUM absolute value (dJ contains the
    % absolute values of the derivatives in the 2 components)
    disp(['The largest shape derivative is ' num2str(max(dJ(:,2)))]); % Consider only the 2nd column, for now
    est_shape_der = [est_shape_der, max(dJ(:,2))];
    if max(dJ(:,2)) < defeaturing_tol % Consider only the 2nd column, for now
        disp('DEFEATURING LOOP: estimator below threshold')
        break;
    end
    
    % Build geometric mesh and space with one level of refinement by the boundary
    [hmsh_geom_ref, hspace_geom_ref] = build_hspace_from_cells(hmsh_geom.ndim, sp_degree, nel_lev0, active_elements_geom_ref, 'active', 'standard', true);
%     hmsh_plot_cells(hmsh_geom_ref, 10, figure);  view(2);
%     title(['Geometric mesh refined boundary ITER ' num2str(iter)],'FontSize',16);
    
    % CHOOSE NEW MESH: select only some of the dJ and build the
    % corresponding mesh: hmsh_new;
    dofs_kept = select_dofs_from_shapeder(dJ, hspace_geom_ref, defeaturing_tol);
    elements_kept = compute_cells_to_refine(hspace_geom_ref, hmsh_geom_ref, dofs_kept); % select elements that should belong to the new geometrical mesh
    el_to_refine = cell(hmsh_geom.nlevels,1);
    for ilev = 1 : hmsh_geom.nlevels % Assumption: hmsh_geom_ref has only one level more than hmsh_geom
        el_to_refine{ilev} = intersect(hmsh_get_parent(hmsh_geom_ref, ilev + 1, elements_kept{ilev + 1}), hmsh_geom.active{ilev});
    end
    % Refine to get the new geometrical mesh
    [hmsh_new, ~] = hmsh_refine (hmsh_geom, el_to_refine);
      
    
    % Update info for the next iteration in python
    active_elements_geom = hmsh_new.active;
%     hmsh_plot_cells(hmsh_new, 10, figure);  view(2);
%     title(['Refined geometric mesh ITER ' num2str(iter)],'FontSize',16);
    
end

disp('DEFEATURING LOOP: Reached maximum number of iterations');

% COMPUTE ERROR ON THE FUNCTIONAL
Jerror_rel = abs(Jexact - Jvalues)./Jexact;
% figure(); plot(ndofs_geom, Jerror_rel, '*-'); grid on; 
% hold on; plot(ndofs_geom, 1./(ndofs_geom.^2), '--')
figure(); semilogy(ndofs_geom, Jerror_rel, '*-'); grid on;
hold on ; semilogy(ndofs_geom, 1000./(ndofs_geom.^2), '--', ndofs_geom, 1./ndofs_geom, '.-')
xlabel('ndofs geometry'); ylabel('rel error $J(\Omega)$','Interpreter','latex');
title('Relative errors on $J(\Omega)$','Fontsize',16,'Interpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define neumann data 
function g = neumann_BC (x, y, ind)
  switch ind
    case 1
      g = 1 + 0.*x.*y;
    case 2
      g = 0.*x.*y;
    case 3
      g = 0.*x.*y;
    case 4
      g = 0.*x.*y;
    otherwise
      error ('g_nmnn: unknown reference number');
  end
end