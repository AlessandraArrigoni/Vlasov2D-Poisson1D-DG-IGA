%% WEAK LANDAU DAMPING TEST with SPLITTING + CRANK NICOLSON TIME SCHEME
% The first subproblem considered is the advection in v-direction, i.e.
% df / dt + E * df / dv = 0
%
% Assumptions: E=grad(\phi), Poisson: \delta\phi = rho_0 - rho with rho_0 background density.
% Polynomials degree for DG method: k; splines order: k+2 with highest
% regularity k+1 on internal and boundary knots (periodic boundary conditions)
%
% The vector nqn stores the number of quadrature nodes per direction: it can be
% different since the polynomial degree is not the same in x and v direction;
% this allows to save some time in assembling the matrices.
%
% The "small" refers to the domain in v (we tried a larger domain [-8,8])

Data= struct(  'testname',         'smallWeakLandau',...
               'domain',           [0,4*pi; -6,6],...     
               'initial_f',        @(x,y) 1/sqrt(2*pi).*(1+0.01*cos(0.5*x)).*exp(-0.5*y.^2),... 
               'f_eq',             @(x,y) 1/sqrt(2*pi).*exp(-0.5*y.^2),... % Maxwellian equilibrium density      
               'rho_0',            1,... % Background density for Poisson
               'fem',              'Q3',...
               'splines',          5,... % Degree to solve Poisson (regularity is deg-1)
               'BC',               'Period',... % Boundary conditions on x-domain
               'nqn',              [6,4],...    % Quadrature nodes in [x,v] direction
               'nref',             5,...   % Refinement level: 2^nref cells per direction
               'time',             0,...
               'dt_scaling',       4,...   % Scaling factor: dt = h / dt_scaling
               'Tend',             100,... % Final simulation time
               'type_mesh',        'CART',...
               'damp_file',        10,...  % Save the solution each n = damp_file "seconds" 
               'computeErrors',    0);  
           
% Call timescheme
filenames = Splitting2orderB(Data);

% Uncomment the following line if run on EPFL server
% exit