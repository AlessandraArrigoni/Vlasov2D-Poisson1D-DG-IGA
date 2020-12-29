%% CONVERGENCE TEST with SPLITTING + CRANK NICOLSON TIME SCHEME
% The first subproblem considered is the advection in x-direction, i.e.
% df / dt + v * df / dx = f
%
% Test in paper by ROSSMANITH SEAL (TIME DEPENDENT)
% Assumptions: E=grad(\phi), Poisson: \delta\phi = rho_0 - rho with rho_0 background density.
% Polynomials degree for DG method: k; splines order: k+2 with highest
% regularity k+1 on internal and boundary knots (periodic boundary conditions)
%
% The vector nqn stores the number of quadrature nodes per direction: it can be
% different since the polynomial degree is not the same in x and v direction;
% this allows to save some time in assembling the matrices.

%% POLYNOMIALS Q1 - SPLINES order 3

Data= struct(  'testname',         'convergence',...
               'domain',           [-pi,pi; -4,4],...              
               'exact_sol',        @(x,y,t) (2-cos(2*x - 2*pi*t)).*exp(-0.25*(4*y-1).^2),...    
               'source',           @(x,y,t) 0.5*exp(-0.25*(4*y-1).^2).*(sin(2*x-2*pi*t)).*((2*sqrt(pi)+1).*(4*y-2*sqrt(pi)) - sqrt(pi).*(4*y-1).*cos(2*x-2*pi*t)),...  
               'rho_0',            sqrt(pi),... % Background density for Poisson
               'fem',              'Q1',...
               'splines',          3,...
               'BC',               'Period',... % Boundary conditions in \Omega_x
               'nqn',              [3,2],...    % Number of quadrature nodes per direction
               'nref',             [3,4,5,6,7],... % Refinement levels (dyadic)
               'time',             0,...
               'dt_scaling',       3,... % Scaling factor: dt = h / dt_scaling
               'Tend',             1,... % Final simulation time
               'type_mesh',        'CART',...
               'computeErrors',    true,...
               'uex_poisson',      @(x,t) cos(2*x-2*pi*t)*sqrt(pi)/8,...
               'gradex_poisson',   @(x,t) -sqrt(pi)/4*sin(2*x-2*pi*t));  
           
% Save a variable with the names of the files storing the results of the
% simulation (one file per refinement level)
filenames = Splitting2orderA_convergence(Data);
save('filenamesAQ1','filenames');

% Uncomment the following line if run on the epfl server
% exit 

%% POLYNOMIALS Q2 - SPLINES order 4

Data= struct(  'testname',         'convergence',...
               'domain',           [-pi,pi; -4,4],...              
               'exact_sol',        @(x,y,t) (2-cos(2*x - 2*pi*t)).*exp(-0.25*(4*y-1).^2),...    
               'source',           @(x,y,t) 0.5*exp(-0.25*(4*y-1).^2).*(sin(2*x-2*pi*t)).*((2*sqrt(pi)+1).*(4*y-2*sqrt(pi)) - sqrt(pi).*(4*y-1).*cos(2*x-2*pi*t)),...  
               'rho_0',            sqrt(pi),... % Background density for Poisson
               'fem',              'Q2',...
               'splines',          4,...
               'BC',               'Period',... 
               'nqn',              [5,3],...    
               'nref',             [3,4,5,6,7],...
               'time',             0,...
               'dt_scaling',       3,... % Scaling factor: dt = h / dt_scaling
               'Tend',             1,...
               'type_mesh',        'CART',...
               'computeErrors',    true,...
               'uex_poisson',      @(x,t) cos(2*x-2*pi*t)*sqrt(pi)/8,...
               'gradex_poisson',   @(x,t) -sqrt(pi)/4*sin(2*x-2*pi*t));  
           
        
filenames = Splitting2orderA_convergence(Data);
save('filenamesAQ2','filenames');

% Uncomment the following line if run on the epfl server
% exit 

%% POLYNOMIALS Q3 - SPLINES order 5

Data= struct(  'testname',         'convergence',...
               'domain',           [-pi,pi; -4,4],...              
               'exact_sol',        @(x,y,t) (2-cos(2*x - 2*pi*t)).*exp(-0.25*(4*y-1).^2),...    
               'source',           @(x,y,t) 0.5*exp(-0.25*(4*y-1).^2).*(sin(2*x-2*pi*t)).*((2*sqrt(pi)+1).*(4*y-2*sqrt(pi)) - sqrt(pi).*(4*y-1).*cos(2*x-2*pi*t)),...  
               'rho_0',            sqrt(pi),... % Background density for Poisson
               'fem',              'Q3',...
               'splines',          5,...
               'BC',               'Period',... 
               'nqn',              [6,4],...    
               'nref',             [3,4,5,6,7],...
               'time',             0,...
               'dt_scaling',       4,... % Scaling factor: dt = h / dt_scaling
               'Tend',             1,...
               'type_mesh',        'CART',...
               'computeErrors',    true,...
               'uex_poisson',      @(x,t) cos(2*x-2*pi*t)*sqrt(pi)/8,...
               'gradex_poisson',   @(x,t) -sqrt(pi)/4*sin(2*x-2*pi*t));  
           
        
filenames = Splitting2orderA_convergence(Data);
save('filenamesAQ3','filenames');

% Uncomment the following line if run on the epfl server
% exit 