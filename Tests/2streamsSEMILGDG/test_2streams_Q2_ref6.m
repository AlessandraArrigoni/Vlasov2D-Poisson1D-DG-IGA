%% TWO STREAMS INSTABILITY TEST with SEMILAGRANGIAN DG SCHEME
% The first subproblem considered is the advection in x-direction, i.e.
% df / dt + v * df / dx = f
%
% Initial condition and parameters as in paper from Ayuso, Hajian (2012)
%
% Assumptions: E=grad(\phi), Poisson: \delta\phi = rho_0 - rho with rho_0 background density.
% Polynomials degree for DG method: k; splines order: k+2 with highest
% regularity k+1 on internal and boundary knots (periodic boundary conditions)
%
% The vector nqn stores the number of quadrature nodes per direction: it can be
% different even if with this method it is not necessary if we choose the same
% polynomial degree for the DG basis in each direction.

Data= struct(  'testname',         'twostreams',...
               'domain',           [0,4*pi; -6,6],...     
               'initial_f',        @(x,y) y.^2./sqrt(8*pi).*(2 - cos(0.5*(x-2*pi))).*exp(-0.5*y.^2),...      
               'rho_0',            1,...        % Background density for Poisson
               'fem',              'Q2',...
               'fem1D',            'P2',...
               'splines',          4,...        % Degree to solve Poisson (regularity is deg-1)
               'BC',               'Period',... % Boundary conditions on x-domain
               'nqn',              [3,3],...    % Quadrature nodes in [x,v] direction
               'nref',             6,...    % Refinement level: 2^nref cells per direction
               'time',             0,...
               'dt_scaling',       4,...    % Scaling factor: dt = h / dt_scaling
               'Tend',             100,...  % Final simulation time
               'type_mesh',        'CART',...
               'damp_file',        15);     % Save the solution each n = damp_file "seconds"  
           
% Call semi-LG method
filenames = semiLGDG(Data);

% Uncomment next line if run on server
% exit