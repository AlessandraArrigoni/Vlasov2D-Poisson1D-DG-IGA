%% ANALYSIS OF THE RESULTS WITH PLOTS AND DECAY RATE ESTIMATE
% Load the file storing the results before running the script

% PLOT DIAGNOSTICS    
tt = linspace(0, data.Tend, length(diagnostics.d_energy_rel));
figure()
subplot(3,1,1)
semilogy(tt,abs(diagnostics.d_energy_rel)); xlabel('t'); ylabel('energy','FontSize',14);
title(['Fem = ' Data.fem ', spline deg = ' num2str(Data.splines)...
    ', cells = ' num2str(2^diagnostics.ref) ],'FontSize',18);
% title(['Conserved properties: fem = ' Data.fem ', spline deg = ' num2str(Data.splines)...
%     ', cells = ' num2str(2^diagnostics.ref) ', \Deltat = \Deltax/' num2str(Data.dt_scaling)],...
%     'FontSize',18);

subplot(3,1,2)
semilogy(tt,abs(diagnostics.d_charge_rel)); xlabel('t'); ylabel('charge','FontSize',14);
ylim([1e-13,1e-9])

subplot(3,1,3)
semilogy(tt,abs(diagnostics.d_L2normVlasov_rel)); xlabel('t'); ylabel('L2norm','FontSize',14);
ylim([1e-15,1e-4]); yticks([1e-15 1e-10 1e-5]);

%% ANALYZE ELECTRIC FIELD DAMPING

% Visualise L2 norm evolution before starting the analysis
tt = linspace(0, Data.Tend, length(diagnostics.L2normE));
[log_max, indices] = findpeaks(log(diagnostics.L2normE));
figure()
semilogy(tt, diagnostics.L2normE, tt(indices), exp(log_max),'v'); xlabel('t'); ylabel('log || E(t) ||_{L^2}');

% Select how many data to consider (might be changed!)
log_max_red = log_max(3:9);
pos_max_red = tt(indices(3:9));
coef = polyfit(pos_max_red, log_max_red,1); % Linear interpolation
damping_factor = coef(1)
damping_teor = -0.1533;

% Frequency: \omega = \pi / time_period since the frequency of the L2 norm is
% twice as large as the frequency of the electric field itself.
% time_period = time_interval / number of peaks contained 
omega_computed = pi*(length(pos_max_red)-1)/(pos_max_red(end)-pos_max_red(1))
omega_teor = 1.4156;

figure()
semilogy(tt, diagnostics.L2normE,'Color','blue','LineWidth',0.8);
hold on
semilogy(tt(indices(3)),exp(log_max(3)),'v','Color','black','MarkerSize',6,'LineWidth',1)
semilogy(tt(indices(9)),exp(log_max(9)),'v','Color','blue','MarkerSize',6,'LineWidth',1)
xlabel('t'); ylabel('log || E(t) ||_{L^2}');
% semilogy(tt(1:floor(length(tt)/2)), exp(coef(1)*tt(1:floor(length(tt)/2))+coef(2)),'-.',...
%     tt(1:floor(length(tt)/2)),exp(damping_teor*tt(1:floor(length(tt)/2))+coef(2)), 'LineWidth',1.3);
% legend('L^2 norm E','Damping computed', 'Damping literature');

% title(['Weak Landau damping: fem = ' Data.fem ', spline deg = ' num2str(Data.splines)...
%     ', cells = ' num2str(2^diagnostics.ref) ', \Deltat = \Deltax /' num2str(Data.dt_scaling)],...
%     'FontSize',18);

semilogy(ttq2(1:floor(length(ttq2)/2)), exp(-0.1533*ttq2(1:floor(length(ttq2)/2))+coef(2)),'r-','LineWidth',1);
legend({'L^2 norm E','First data point','Last data point',['\gamma = ' num2str(damping_factor)]},'FontSize',12,'Location','SE' );

title(['Fem = ' Data.fem ', spline deg = ' num2str(Data.splines) ', cells = ' num2str(2^diagnostics.ref)],...
    'FontSize',18);

%% SAVE RESULTS TO VTK

filename = ['landau_semiLGDG_' Data.fem '_ref' num2str(ref) '_t'];
tempname = [filename '0'];
[basis_plot, coord_plot] = write_solDG_to_vtk(tempname, Data, femregion, ff.first, [4,8],1);

for j = 1:10
    tempname = [filename num2str(j)];
    tempsol = ff.(['n' num2str(j)]);
    [basis_plot, coord_plot] = write_solDG_to_vtk(tempname,Data, femregion, tempsol, [4,8],1, basis_plot, coord_plot);
end

fprintf('\nEND saving data for paraview\n');

%% Plots used in the report and slides 
% Parameters: 
% plotboxaspect ratio = 1 , 0.8025, 0.8025,
% outer position = 0.01042027729636,0,0.498694306036973,1
% position = 0.0753,0.11,0.423,0.815

figure() 
markers{6} = '-'; markers{7} = '-.';
colors{6} = 'blue'; colors{7} = 'black';
for ref = 6:7
    load(['smallWeakLandau_semiLGDG_Q3_ref' num2str(ref) '.mat']);
    tt = linspace(0, Data.Tend, length(diagnostics.d_energy_rel));

    subplot(3,2,[1,3,5]) % Landau damping
    semilogy(tt, diagnostics.L2normE ,markers{ref},'color',colors{ref}, 'LineWidth',1);
    hold on
    title('\textbf{Polynomial degree r = 3 with splines degree 4}','Interpreter','latex','FontSize',20);
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('$\log || E(t) ||_{L^2}$','Interpreter','latex','FontSize',16);
    
    
    subplot(3,2,2)
    hold on
    semilogy(tt,abs(diagnostics.d_energy_rel),markers{ref},'color',colors{ref},'Linewidth',1); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('Energy','Interpreter','latex','FontSize',16);
    ylim([-1e-7, 5*1e-7]);
    
    subplot(3,2,4)
    hold on
    semilogy(tt,abs(diagnostics.d_charge_rel),markers{ref},'color',colors{ref},'Linewidth',0.8); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('Charge','Interpreter','latex','FontSize',16);
    ylim([0, 2*1e-10]);
    
    subplot(3,2,6)
    hold on
    semilogy(tt,abs(diagnostics.d_L2normVlasov_rel),markers{ref},'color',colors{ref},'Linewidth',1); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('$L^2$ norm','Interpreter','latex','FontSize',16);
    ylim([-1e-5, 3*1e-5]);
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

%% Plot again for the slides.

figure() 
markers{6} = '-'; markers{7} = '-.';
colors{6} = 'blue'; colors{7} = 'black';
for ref = 6:7
    load(['smallWeakLandau_semiLGDG_Q3_ref' num2str(ref) '.mat']);
    tt = linspace(0, Data.Tend, length(diagnostics.d_energy_rel));

    subplot(3,1,1)
    hold on
    semilogy(tt,abs(diagnostics.d_energy_rel),markers{ref},'color',colors{ref},'Linewidth',1); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('Energy','Interpreter','latex','FontSize',16);
    ylim([-1e-7, 5*1e-7]);
    
    subplot(3,1,2)
    hold on
    semilogy(tt,abs(diagnostics.d_charge_rel),markers{ref},'color',colors{ref},'Linewidth',0.8); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('Charge','Interpreter','latex','FontSize',16);
    ylim([0, 2*1e-10]);
    
    subplot(3,1,3)
    hold on
    semilogy(tt,abs(diagnostics.d_L2normVlasov_rel),markers{ref},'color',colors{ref},'Linewidth',1); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('$L^2$ norm','Interpreter','latex','FontSize',16);
    ylim([-1e-5, 3*1e-5]);
end