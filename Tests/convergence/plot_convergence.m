%% CONVERGENCE PLOTS
% Error vs meshsize at final time for the test with Runge Kutta time scheme

errL2_vlasovq1 = []; errL2_vlasovq2 = []; errL2_vlasovq3 = [];
errL2_poissonq1 = []; errL2_poissonq2 = []; errL2_poissonq3 = [];
errL2_electricq1 = []; errL2_electricq2 = []; errL2_electricq3 = [];
hh = [];

for r = 3:7
    load(['convergence_rk4_Q1_ref' num2str(r) '.mat']);
    errL2_vlasovq1 = [errL2_vlasovq1, vlasov_L2];
    errL2_poissonq1 = [errL2_poissonq1, poisson_l2];
    errL2_electricq1 = [errL2_electricq1, electric_l2];
    hh = [hh, femregion.h];
end
for r = 3:7
    load(['convergence_rk4_Q2_ref' num2str(r) '.mat']);
    errL2_vlasovq2 = [errL2_vlasovq2, vlasov_L2];
    errL2_poissonq2 = [errL2_poissonq2, poisson_l2];
    errL2_electricq2 = [errL2_electricq2, electric_l2];
end
for r = 3:7
    load(['convergence_rk4_Q3_ref' num2str(r) '.mat']);
    errL2_vlasovq3 = [errL2_vlasovq3, vlasov_L2];
    errL2_poissonq3 = [errL2_poissonq3, poisson_l2];
    errL2_electricq3 = [errL2_electricq3, electric_l2];
end

% Compute relative errors
load('exactL2norms.mat'); % Load the norms of the exact solution computed with many quadrature nodes
errL2_vlasovq1 = errL2_vlasovq1/L2normexactVlasov; errL2_vlasovq2 = errL2_vlasovq2/L2normexactVlasov; errL2_vlasovq3 = errL2_vlasovq3/L2normexactVlasov;
errL2_poissonq1 = errL2_poissonq1/L2uex; errL2_poissonq2 = errL2_poissonq2/L2uex; errL2_poissonq3 = errL2_poissonq3/L2uex; 
errL2_electricq1 = errL2_electricq1/L2graduex; errL2_electricq2 = errL2_electricq2/L2graduex; errL2_electricq3 = errL2_electricq3/L2graduex; 

%% start plotting 
figure()
subplot(1,3,1)
loglog(hh, errL2_vlasovq1, 'v-','LineWidth',1,'MarkerSize',6)
hold on; grid on
loglog(hh, errL2_vlasovq2, '*-', 'LineWidth',1,'MarkerSize',6)
loglog(hh, errL2_vlasovq3, 'o-', 'LineWidth',1,'MarkerSize',6)
loglog( hh, 100*hh.^2, '-.','LineWidth',1,'Color','k')
loglog(hh, 200*hh.^3 ,'-','LineWidth',1,'Color','k')
loglog(hh, 200*hh.^4, '--','LineWidth',1,'Color','k')
legend({'r = 1','r = 2','r = 3','$O(h^2)$','$O(h^3)$','$O(h^4)$'},'Interpreter','latex', 'FontSize',12,'Location','SE')
xlabel('h','FontSize',12); ylabel('$\| f- f_h\|_{L^2} \, / \, \|f\|_{L^2}$','Interpreter','latex','FontSize',12)


subplot(1,3,2)
loglog(hh, errL2_poissonq1, 'v-','LineWidth',1,'MarkerSize',6)
hold on; grid on
loglog(hh, errL2_poissonq2, '*-', 'LineWidth',1,'MarkerSize',6) 
loglog(hh, errL2_poissonq3, 'o-',  'LineWidth',1,'MarkerSize',6)
loglog( hh, 10*hh.^2,'-.','LineWidth',1,'Color','k');
loglog(hh, 5*hh.^4 ,'--','LineWidth',1,'Color','k');
legend({'splines deg 3','splines deg 4','splines deg 5','$O(h^2)$','$O(h^4)$'},'Interpreter','latex', 'FontSize',12,'Location','SE')
xlabel('h','FontSize',12); ylabel('$\| \phi- \phi_h\|_{L^2} \, / \, \|\phi\|_{L^2}$','Interpreter','latex','FontSize',12)

subplot(1,3,3)
loglog(hh, errL2_electricq1, 'v-', 'LineWidth',1,'MarkerSize',6)
hold on; grid on
loglog(hh, errL2_electricq2, '*-',  'LineWidth',1,'MarkerSize',6) 
loglog(hh, errL2_electricq3, 'o-',  'LineWidth',1,'MarkerSize',6)
loglog(hh, 10*hh.^2,'-.','LineWidth',1, 'Color','k')
loglog(hh, 5*hh.^4 ,'--','LineWidth',1,'Color','k')
legend({'splines deg 3','splines deg 4','splines deg 5','$O(h^2)$','$O(h^4)$'},'Interpreter','latex', 'FontSize',12,'Location','SE')
xlabel('h','FontSize',12); ylabel('$\| E- E_h\|_{L^2} \, / \, \|E\|_{L^2}$','Interpreter','latex','FontSize',12)