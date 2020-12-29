%% TEST L2 PROJECTION 1D with SEMILAGRANGIAN DG METHOD
% Script to test the implementation of the MASS MATRIX for the DG method 1D:
% it should be block diagonal and we check the convergence of the L2 projection
% of a regular function on thE discontinuous spaces of polynomials.

data = struct('domain',      [0,1 ; -4,4],...
              'initial_f',   @(x,y) cos(2*pi*x),... % period = pi
              'trasp',       1,...      % Transport coefficient
              'fem',         'Q3',...
              'fem1D',       'P3',...
              'BC',          'Period',...
              'nqn',         [4,4],...  % Enough to build matrices for P3 (integrate polynomials of degree 6)
              'time',         0,...
              'Tend',         0.5,...   % Final simulation time
              'type_mesh',    'CART');
        
nref = [3, 4, 5, 6];
errL2sys = [];
errL2inv = [];

for ref = nref 
    % Mesh and DG finite element space data structures
    [region] = generate_mesh(data, ref); % region 2D
    [femregion] = create_dof(data, region); % FE space 2D
    femregionX = create_dof1D(data, region, 'x'); % FE space 1D in x
    [nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(data.nqn);
    
    % Evaluate basis 1D
    shapeX = basis_lagrange1D(data.fem1D);
    [phiX, gradX] = evalshape1D(shapeX, nodes_1D{1}, femregionX.nln); 
    basisX = struct('nodes_1D', nodes_1D{1}, 'w_1D', w_1D{1}, ...
            'dphiq', phiX, 'gradphi', gradX);
    
    % Mass matrix 1D
    [M, Minv] = mass_matrix1D(femregionX, basisX);
    f = rhs_givenf1D(femregionX, @(x) data.initial_f(x,0) , basisX);

    % Solution of the linear system
    uh1 = M\f;
    uh2 = Minv*f;

    % Post-processing
    title1 = ['L2 projection (solving the system) with fem ' data.fem1D ' ref ' num2str(ref)];
    postprocessing1D(femregionX, data, uh1, 0, 6, title1);
    title2 = ['L2 projection (computing the inverse) with fem ' data.fem1D ' ref ' num2str(ref)];
    postprocessing1D(femregionX, data, uh2, 0, 6, title2);
    
    % Errors
    [err1]= compute_errors1D(data, femregionX, uh1, 0);
    errL2sys = [errL2sys, err1.E_L2];
    [err2] = compute_errors1D(data, femregionX, uh2, 0);
    errL2inv = [errL2inv, err2.E_L2];
end

% Analysis
hh = 0.5.^nref;
figure()
loglog(hh, errL2sys, '*-', hh, hh.^(femregion.degree) , hh, hh.^(femregion.degree + 1))
hold on
loglog(hh, errL2inv, 'o--','MarkerSize',6)
legend('Error L2 sys', ['h^' num2str(femregion.degree)], ['h^' num2str(femregion.degree + 1)], 'Error L2 inv')
title(['Error con fem = ' data.fem])