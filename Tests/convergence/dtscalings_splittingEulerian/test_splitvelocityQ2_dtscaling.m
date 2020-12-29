%% Splitting II order starting from transport in x (advection = velocity)

% Assumptions: E=grad(\phi), 
% Poisson: \delta\phi = rho_0 - rho with rho_0 background density.
% Polynomials degree for DG method: k; splines order: k+2 with highest
% regularity k+1 on internal and boundary knots (periodic boundary conditions)
Data= struct(  'testname',         'conv_dtscaling',...
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
               'scalings',         [2,4,8,16],... % scaling factor: dt = h / dt_scaling
               'Tend',             1,...
               'type_mesh',        'CART',...
               'computeErrors',    true,...
               'uex_poisson',      @(x,t) cos(2*x-2*pi*t)*sqrt(pi)/8,...
               'gradex_poisson',   @(x,t) -sqrt(pi)/4*sin(2*x-2*pi*t));  
           
        
filenames = Splitting2orderA_conv_scaling(Data);

save('filenamesAQ2','filenames');

% Uncomment the following line if run on server
% exit 


%% Analyze errors
% CFL is defined as max(trasp)dt/dx2D = max(trasp)/dt_scaling where
% max(trasp) = 4 given the domain we consider

% For each refinement we saved vectors containing the errors for each dt_scaling
nref = [3,4,5,6,7];
scalings = [2,4,8,16];
vlasov = zeros(length(nref), length(scalings));
poisson = zeros(length(nref), length(scalings));
electric = zeros(length(nref), length(scalings));
hh = zeros(length(nref),1);

for i = 1:length(nref)
    ref = nref(i);
    load(filenames.(['ref' num2str(ref)]));
    
    vlasov(i,:) = errVlasov;
    poisson(i,:) = errPoisson;
    electric(i,:) = errElectric;    
    hh(i) = femregion.h;
end
% Compute relative errors
load('exactL2norms.mat');
vlasov = vlasov / L2normexactVlasov;
poisson = poisson / L2uex;
electric = electric / L2graduex;
%% Plot errors
colors = hsv(length(scalings)); % prism hsv jet
markers = ['o','+','*','d']; % knowing that we used 4 different dt_scalings
figure()


for i = 1:length(Data.scalings)
    subplot(1,3,1)
    loglog(hh, vlasov(:,i),'-','Marker',markers(i),'Color',colors(i,:),'LineWidth',1,'MarkerSize',7)
    hold on; grid on 
    leg1{i} = ['CFL = ' num2str(4/Data.scalings(i),'%.2f')];
    
    subplot(1,3,2)
    loglog(hh, poisson(:,i),'-','Marker',markers(i),'Color',colors(i,:), 'LineWidth',1,'MarkerSize',7)
    hold on; grid on 
    leg2{i} = ['CFL = ' num2str(4/Data.scalings(i),'%.2f')];
    
    subplot(1,3,3)
    loglog(hh, electric(:,i),'-','Marker',markers(i),'Color',colors(i,:), 'LineWidth',1,'MarkerSize',7)
    hold on; grid on 
    leg3{i} = ['CFL = ' num2str(4/Data.scalings(i),'%.2f')];
end

% Q2
subplot(1,3,1); 
loglog(hh, 150*hh.^(femregion.degree + 1),'k--','LineWidth',1); 
leg1{i+1} = ['$O(h^' num2str(femregion.degree+1) ')$'];
legend(leg1,'Interpreter','latex','Location','SE','FontSize',14); xlabel('h','FontSize',14); 
ylabel('$\| f- f_h\|_{L^2} \, / \, \|f\|_{L^2}$','Interpreter','latex','FontSize',12);

subplot(1,3,2); 
loglog(hh, hh.^(femregion.degree),'k-.','LineWidth',1); 
loglog(hh, 5*hh.^(femregion.degree+2),'k--','LineWidth',1); 
leg2{i+1} = ['$O(h^' num2str(femregion.degree) ')$'];
leg2{i+2} = ['$O(h^' num2str(femregion.degree+2) ')$'];
legend(leg2,'Interpreter','latex','Location','SE','FontSize',14); xlabel('h','FontSize',14); 
ylabel('$\| \phi- \phi_h\|_{L^2} \, / \, \|\phi\|_{L^2}$','Interpreter','latex','FontSize',12)
title(['\textbf{Polynomial degree r = ' num2str(femregion.degree) ' with splines degree ' num2str(femregion.degree+2) '}'],'Interpreter','latex','FontSize',20) 


subplot(1,3,3);
loglog(hh, 0.5*hh.^(femregion.degree),'k-.','LineWidth',1); 
loglog(hh, 10*hh.^(femregion.degree+2),'k--','LineWidth',1); 
leg3{i+1} = ['$O(h^' num2str(femregion.degree) ')$'];
leg3{i+2} = ['$O(h^' num2str(femregion.degree+2) ')$'];
legend(leg3, 'Interpreter','latex','Location','SE','FontSize',14); xlabel('h','FontSize',14); 
ylabel('$\| E- E_h\|_{L^2} \, / \, \|E\|_{L^2}$','Interpreter','latex','FontSize',12)
