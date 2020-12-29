%% ANALYSIS OF THE RESULTS WITH PLOTS 
% Load the file storing the results before running the script

% PLOT DIAGNOSTICS
tt = linspace(0, Data.Tend, length(diagnostics.d_energy_rel));
figure()
subplot(3,1,1)
semilogy(tt,abs(diagnostics.d_energy_rel)); xlabel('t'); ylabel('energy','FontSize',14);
yticks([1e-10,1e-5])
title(['Fem = ' Data.fem ', spline deg = ' num2str(Data.splines)...
    ', cells = ' num2str(2^diagnostics.ref) ],'FontSize',18);
% title(['Conserved properties: fem = ' Data.fem ', spline deg = ' num2str(Data.splines)...
%     ', cells = ' num2str(2^diagnostics.ref) ', \Deltat = \Deltax/' num2str(Data.dt_scaling)],...
%     'FontSize',18);

subplot(3,1,2)
semilogy(tt,abs(diagnostics.d_charge_rel)); xlabel('t'); ylabel('charge','FontSize',14);
ylim([1e-16, 1e-11]); yticks([1e-15,1e-13, 1e-11]);

subplot(3,1,3)
semilogy(tt,abs(diagnostics.d_L2normVlasov_rel)); xlabel('t'); ylabel('L2norm','FontSize',14);
yticks([1e-10,1e-6,1e-2]);


%% SAVE RESULTS TO VTK
filename = ['2streams_splittingA_' Data.fem '_ref' num2str(ref) '_t'];
tempname = [filename '0'];
[basis_plot, coord_plot] = write_solDG_to_vtk(tempname, Data, femregion, ff.first, [7,7],0);

for j = 1:6 % NB: check they are in fact 6
    tempname = [filename num2str(j)];
    tempsol = ff.(['n' num2str(j)]);
    [basis_plot, coord_plot] = write_solDG_to_vtk(tempname,Data, femregion, tempsol, [7,7],0, basis_plot, coord_plot);
end

tempname = [filename '7'];
[basis_plot, coord_plot] = write_solDG_to_vtk(tempname,Data, femregion, ff.last, [7,7],0, basis_plot, coord_plot);
fprintf('\nEND saving data for paraview\n');

%% Plot diagnostics for latex 
% Q2
figure() 
markers{6} = '-'; markers{7} = '-.';
colors{6} = 'blue'; colors{7} = 'black';
for ref = 6:7
    load(['twostreams_splittingA_Q2_ref' num2str(ref) '.mat']);
    tt = linspace(0, Data.Tend, length(diagnostics.d_energy_rel));
    
    subplot(3,1,1)
    hold on
    plot(tt,abs(diagnostics.d_energy_rel),markers{ref},'color',colors{ref},'Linewidth',1); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('Energy','Interpreter','latex','FontSize',16);
    ylim([-1e-5, 4*1e-4]);
    title('\textbf{Polynomial degree r = 2 with splines degree 4}','Interpreter','latex','FontSize',20);
    
    
    subplot(3,1,2)
    hold on
    plot(tt,abs(diagnostics.d_charge_rel),markers{ref},'color',colors{ref},'Linewidth',0.8); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('Charge','Interpreter','latex','FontSize',16);
    ylim([-1e-14, 5*1e-13]);
    
    subplot(3,1,3)
    hold on
    plot(tt,abs(diagnostics.d_L2normVlasov_rel),markers{ref},'color',colors{ref},'Linewidth',1); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('$L^2$ norm','Interpreter','latex','FontSize',16);
    ylim([0,0.15]);
    
end
legend({'64 cells','128 cells'},'Interpreter','latex','FontSize',14)

%% Plot diagnostics for latex
% Q3
figure() 
markers{6} = '-'; markers{7} = '-.';
colors{6} = 'blue'; colors{7} = 'black';
for ref = 6:7
    load(['twostreams_splittingA_Q3_ref' num2str(ref) '.mat']);
    tt = linspace(0, Data.Tend, length(diagnostics.d_energy_rel));
    
    subplot(3,1,1)
    hold on
    plot(tt,abs(diagnostics.d_energy_rel),markers{ref},'color',colors{ref},'Linewidth',1); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('Energy','Interpreter','latex','FontSize',16);
    ylim([-1e-5, 4*1e-4]);
    title('\textbf{Polynomial degree r = 3 with splines degree 5}','Interpreter','latex','FontSize',20);
    
    
    subplot(3,1,2)
    hold on
    plot(tt,abs(diagnostics.d_charge_rel),markers{ref},'color',colors{ref},'Linewidth',0.8); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('Charge','Interpreter','latex','FontSize',16);
    ylim([-1e-14, 5*1e-13]);
    
    subplot(3,1,3)
    hold on
    plot(tt,abs(diagnostics.d_L2normVlasov_rel),markers{ref},'color',colors{ref},'Linewidth',1); 
    xlabel('t','Interpreter','latex','FontSize',16); ylabel('$L^2$ norm','Interpreter','latex','FontSize',16);
    ylim([0,0.15]);
    
end
legend({'64 cells','128 cells'},'Interpreter','latex','FontSize',14)