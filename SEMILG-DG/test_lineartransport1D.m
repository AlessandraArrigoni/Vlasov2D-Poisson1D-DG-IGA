%% TEST LINEAR TRANSPORT EQUATION 1D with SEMILAGRANGIAN DG METHOD
% Script to test the implementation of the SEMILG-DG method on the linear 
% transport equation 1D with periodic or "compact support" (i.e. the solution
% is 0 outside the domain) boundary conditions.
% 
% The problem we consider is df/dt + trasp df/dx = 0
% To understand the meaning of "CIRCULANT MATRIX" refer to the paper by
% Crouseilles et al. (2011)
%
% We compute the errors in L2 norm with respect to the exact solution 
% (translation of the initial profile)
% NB: the timestep must be sufficiently small to ensure that the
% characteristics' starting point falls in the neighbouring cell at most.
% 

% TESTS : sol = sin(2*pi*x) on (0,1) with periodic BC or
% sol = exp(-100*(x-1).^2) on (0,2) with "compact support" BC 

data = struct('domain',      [0,1;-1,1],...
              'initial_f',   @(x,y) sin(2*pi*x),... % period = pi
              'trasp',       1,...     % Transport coefficient in the equation
              'fem',         'Q1',...
              'fem1D',       'P1',...
              'BC',          'Period',...
              'nqn',         [2,2],... % Number of quadrature nodes in (x,v)
              'time',         0,...
              'Tend',         0.5,...  % Final simulation time
              'type_mesh',    'CART');
          
nref = [3,4,5,6,7]; % Refinement levels: 2^nref cells per direction
dt_scaling = 5;     % It must be dt < dx/trasp;
errL1 = [];
errL2 = [];
hh = [];
[nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(data.nqn);
% Translate nodes on the reference element (0,1)
nodes_ref = (nodes_1D{1} + 1)*0.5;

for ref = nref % Loop over refinements
    region = generate_mesh(data, ref);      % Region 2D
    femregion = create_dof(data, region);   % FE space 2D
    femregionX = create_dof1D(data, region, 'x'); % FE space 1D in x
    femregionY = create_dof1D(data, region, 'y'); % FE space 1D in y
    
    dt = femregionX.h/dt_scaling;
    
    % Evaluate basis 1D
    shapeX = basis_lagrange1D(data.fem1D);
    [phiX, gradX] = evalshape1D(shapeX, nodes_ref, femregionX.nln); 
    basisX = struct('nodes_1D', nodes_1D{1}, 'w_1D', w_1D{1}, ...
                'dphiq', phiX, 'gradphi', gradX);

    % Mass matrix 1D
    [M, Minv, Mloc] = mass_matrix1D(femregionX, basisX);

    % Circulant matrix (transport)
    if femregionX.direction == 'x' % Periodic boundary conditions
        [rows, cols] = circulant_matrix1D_structure(femregionX, data.trasp);
    elseif femregionX.direction == 'y' % "Compact support" boundary conditions
        [rows, cols] = circulant_matrix1D_BC_structure(femregionX, data.trasp);
    end
    values = transport_matrix1D(femregionX, shapeX, nodes_ref, w_1D{1}, data.trasp, dt);
    L = sparse(rows, cols, values);

    % Initialize the system
    fold = zeros(femregionX.ndof,1);
    fnew = fold;
    fold = data.initial_f(femregionX.dof, 0);
    postprocessing1D(femregionX, data, fold, 0, 6 ,'Initial condition'); % Plot the initial condition
    
    % Time loop
    for t = 0 : dt : data.Tend-dt
        fnew = Minv *(L * fold);

%         mytitle = ['Solution at t = ' num2str(t+dt)];
%         npoints_plot_4elem = 16;
%         postprocessing1D(femregionX, dati, fnew, t+dt, npoints_plot_4elem, mytitle);
%         pause(0.1)

        fold = fnew;
    end
    postprocessing1D(femregionX, data, fnew, t+dt, 16, 'Solution at t = 0.5');
    errors = compute_errors1D(data, femregionX, fnew, t+dt);
    errL1 = [errL1, errors.E_L1];
    errL2 = [errL2, errors.E_L2];
    hh = [hh, femregionX.h];
    
end % Loop over refinements

%% Analyse convergence

figure()
loglog(hh, errL1, '+-', hh, errL2, 'o-', hh, hh.^(femregionX.degree + 1),'--')
legend('L1 norm','L2 norm',['h^' num2str(femregionX.degree + 1)]);
title(['Errors with polynomials ' femregionX.fem])

pL1 = log(errL1(2:end)./errL1(1:end-1))./log(hh(2:end)./hh(1:end-1))
pL2 = log(errL2(2:end)./errL2(1:end-1))./log(hh(2:end)./hh(1:end-1))
