function [J, filename_out] = exact_geometry_and_functional(filename, exactgeom_thresh, nel_lev0, sp_degree, proj_type, sin_parameters, ...
                                                problem_data, adaptivity_data, Jstring, Jdomain, hmsh_base, hspace_base, plot_flag )
                                            
% EXACT_GEOMETRY_AND_FUNCTIONAL : approximate the geometry with high
% accuracy (exactgeom_thresh) and compute the solution to the pde with a
% very refined mesh to obtain the "exact" value of the considered
% functional.
%
% OUTPUT: 
%          J: "exact" value of the functional on the "exact" geometry
%          filename_out: file that stores info for the "exact" geometry
%
%

if nargin < 13
    plot_flag = false;
end
    
% Check if the file with the exact geometry exists and load it
filename_out = [filename(1:end-4) '_out.mat'];
if ~exist(filename_out,'file')
    active_elements_geom = {[1:prod(nel_lev0)]'};
    save(filename, 'active_elements_geom','nel_lev0','sp_degree','proj_type','sin_parameters','exactgeom_thresh');
    disp([newline 'COMPUTING EXACT GEOMETRY']);
    py.inout.exact_geometry(filename) % call python function
end
load(filename_out);
if ~iscell(active_elements_geom_nofold)
    active_elements_geom_nofold = {active_elements_geom_nofold};
end 

% CONSTRUCT MESH AND SPACE AND COMPUTE THE EXACT SOLUTION (ADAPTIVE)
nlevels = numel(active_elements_geom_nofold)-1; 
adaptivity_data_geom.flag = 'elements';
hmsh_geom = hmsh_base; hspace_geom = hspace_base;
for ref = 1:nlevels
    marked_geom = cell (ref,1);
    marked_geom{ref} = setdiff(hmsh_geom.active{hmsh_geom.nlevels}, active_elements_geom_nofold{ref}); 
    [hmsh_geom, hspace_geom] = adaptivity_refine (hmsh_geom, hspace_geom, marked_geom, adaptivity_data_geom); 
end

% Geometry reconstruction
C = hspace_subdivision_matrix (hspace_geom);
fine_cpoints =  C{hspace_geom.nlevels} * cpoints; 
fine_cpoints = reshape(fine_cpoints', [hmsh_geom.rdim, hspace_geom.space_of_level(end).ndof_dir]);
% Build the geometry mappings
fine_knots = hspace_geom.space_of_level(end).knots;
nurbs = nrbmak(fine_cpoints, fine_knots); 
geometry = geo_load (nurbs);
% Update the mapping
hmsh = hmsh_geom; hspace = hspace_geom;
hmsh = update_mappings(hmsh, geometry);
if plot_flag
    hmsh_plot_cells(hmsh,10,figure); view(2); title('Exact geometry','FontSize',16);
end

% Adaptive loop
it_adapt = 0; gest = [];
while true
    it_adapt = it_adapt + 1;    
    if (~hspace_check_partition_of_unity (hspace, hmsh))
        disp('ERROR: The partition-of-the-unity property does not hold.')
        break
    end

    % Solve original PDE
    u = adaptivity_solve_laplace (hmsh, hspace, problem_data);
    nel(it_adapt) = hmsh.nel; ndof(it_adapt) = hspace.ndof;
    % Estimate and check
    est = adaptivity_estimate_laplace (u, hmsh, hspace, problem_data, adaptivity_data);
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

if plot_flag
    figure()
    sp_plot_solution(u, hspace, geometry, 201); view(2); colorbar();
    title('Exact solution', 'FontSize',16);
end

% Compute the functional
J = compute_functionalJ(Jstring, Jdomain, u, hmsh, hspace);


end

