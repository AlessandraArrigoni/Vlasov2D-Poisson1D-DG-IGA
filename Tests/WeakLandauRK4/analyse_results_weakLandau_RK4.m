%% ANALYSIS OF THE RESULTS WITH PLOTS AND DECAY RATE ESTIMATE
% Load the file storing the results before running the script

%% Plot to compare small and large domain 
% Input files with large domain are not provided but can be easily obtained by
% modifying the existing ones
load('smallWeakLandau_rk4_Q2_ref6.mat')
ttq2 = linspace(0, Dati.Tend, length(diagnostics.L2normE));
smallq2 = diagnostics.L2normE;
[log_maxq2, indicesq2] = findpeaks(log(diagnostics.L2normE));
log_max_red = log_maxq2(3:9);
pos_max_red = ttq2(indicesq2(3:9));
coefq2 = polyfit(pos_max_red, log_max_red,1); % Linear interpolation
load('largeWeakLandau_rk4_Q2_ref6.mat')
largeq2 = diagnostics.L2normE;

% case Q2
figure() 
semilogy(ttq2, smallq2,'Color','blue','LineWidth',1);
hold on
semilogy(ttq2, largeq2, '-.','Color' ,'black','LineWidth',1);
semilogy(ttq2(indicesq2(3)),exp(log_maxq2(3)),'v','Color','red','MarkerSize',8,'LineWidth',1)
semilogy(ttq2(indicesq2(9)),exp(log_maxq2(9)),'^','Color','red','MarkerSize',8,'LineWidth',1)
semilogy(ttq2(1:floor(length(ttq2)/2)), exp(-0.1533*ttq2(1:floor(length(ttq2)/2))+coefq2(2)),'r-','LineWidth',1);
legend({'$\Omega_v$ small', '$\Omega_v$ large','third peak','ninth peak','$\gamma = 0.1533$'},'Interpreter','latex','FontSize',14)
title('\textbf{Polynomial degree r = 2, 64 cells}','Interpreter','latex','FontSize',20)
ylim([1e-7,0.2])
xlabel('t','FontSize',16); ylabel('$\log || E(t) ||_{L^2}$','Interpreter','latex','FontSize',16);

% case Q3
load('smallWeakLandau_rk4_Q3_ref6.mat')
ttq3 = linspace(0, Dati.Tend, length(diagnostics.L2normE));
smallq3 = diagnostics.L2normE;
[log_maxq3, indicesq3] = findpeaks(log(diagnostics.L2normE));
log_max_red = log_maxq2(3:9);
pos_max_red = ttq2(indicesq2(3:9));
coefq3 = polyfit(pos_max_red, log_max_red,1); % Linear interpolation
load('largeWeakLandau_rk4_Q3_ref6.mat')
largeq3 = diagnostics.L2normE;

figure() 
semilogy(ttq3, smallq3,'Color','blue','LineWidth',1);
hold on
semilogy(ttq3, largeq3, '-.','Color' ,'black','LineWidth',1);
semilogy(ttq3(indicesq3(3)),exp(log_maxq3(3)),'v','Color','red','MarkerSize',8,'LineWidth',1)
semilogy(ttq3(indicesq3(9)),exp(log_maxq3(9)),'^','Color','red','MarkerSize',8,'LineWidth',1)
semilogy(ttq3(1:floor(length(ttq3)/2)), exp(-0.1533*ttq3(1:floor(length(ttq3)/2))+coefq3(2)),'r-','LineWidth',1);
legend({'$\Omega_v$ small', '$\Omega_v$ large','third peak','ninth peak','$\gamma = 0.1533$'},'Interpreter','latex','FontSize',14)
title('\textbf{Polynomial degree r = 3, 64 cells}','Interpreter','latex','FontSize',20)
ylim([1e-7,0.2])
xlabel('t','FontSize',16); ylabel('$\log || E(t) ||_{L^2}$','Interpreter','latex','FontSize',16);

%% Plot to compare different mesh size for a given polynomial degree
% Q2 (do the same for Q3, with the corresponding parameters)
figure(1) 
markers{6} = '-'; markers{7} = '-.';
colors{6} = 'blue'; colors{7} = 'black';
for ref = 6:7
    load(['smallWeakLandau_rk4_Q2_ref' num2str(ref) '.mat']);
    tt = linspace(0, Dati.Tend, length(diagnostics.d_energy_rel));

    subplot(3,2,[1,3,5]) % landau damping
    semilogy(tt, diagnostics.L2normE ,markers{ref},'color',colors{ref}, 'LineWidth',1);
    hold on
    title('\textbf{Polynomial degree r = 2 with splines degree 4}','Interpreter','latex','FontSize',20);
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('$\log || E(t) ||_{L^2}$','Interpreter','latex','FontSize',16);
    
    
    subplot(3,2,2)
    hold on
    plot(tt,abs(diagnostics.d_energy_rel),markers{ref},'color',colors{ref},'Linewidth',1); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('Energy','Interpreter','latex','FontSize',16);
    ylim([-1e-14, 5*1e-13]);
    
    subplot(3,2,4)
    hold on
    plot(tt,abs(diagnostics.d_charge_rel),markers{ref},'color',colors{ref},'Linewidth',0.8); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('Charge','Interpreter','latex','FontSize',16);
    ylim([-1e-14, 1e-13]);
    
    subplot(3,2,6)
    hold on
    plot(tt,abs(diagnostics.d_L2normVlasov_rel), markers{ref},'color',colors{ref},'Linewidth',1); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('$L^2$ norm','Interpreter','latex','FontSize',16);
    ylim([-1e-6, 5*1e-6]);
end

subplot(3,2,[1,3,5])
[log_max, indices] = findpeaks(log(diagnostics.L2normE));
log_max_red = log_max(3:9);
pos_max_red = tt(indices(3:9));
coef = polyfit(pos_max_red, log_max_red,1); % Linear interpolation
semilogy(tt(indices(3)),exp(log_max(3)),'v','Color','red','MarkerSize',8,'LineWidth',1.2)
semilogy(tt(indices(9)),exp(log_max(9)),'^','Color','red','MarkerSize',8,'LineWidth',1.2)
semilogy(tt(1:floor(length(tt)/2)), exp(-0.1533*tt(1:floor(length(tt)/2))+0.95*coef(2)),'r-','LineWidth',1);
legend({'64 cells', '128 cells','third peak','ninth peak','$\gamma = 0.1533$'},'Interpreter','latex','FontSize',14)

