%% Convergence check (NOT USED IN THE VLASOV-POISSON PROJECT)
function [pL2, pH1] = convergence_check_IGA(subs, errL2, errH1, degree, contBC)
% Computes the convergence order p and shows the errors convergence plots.
% The input contBC (<degree) can be a scalar or a vector if we want to compare different
% values for the same spline degree. In the second case, errL2 and errH1 are matrices
% with nrows = length(contBC)

h = 1./subs;
pL2 = [];
pH1 = [];
figure()
for i = 1:length(contBC)
    
    subplot(2,1,1)
    loglog(h, errL2(i,:), '*-', h, h.^(degree+1))
    legend('L2 error', ['h^' num2str(degree+1)])
    grid on
    title(['Degree = ' num2str(degree) ' and continuityBC = ' num2str(contBC(i))])
    subplot(2,1,2)
    loglog(h, errH1(i,:), '*-', h, h.^degree)
    legend('H1 error',['h^' num2str(degree)])
    grid on
    
    
    % save convergence rates: matrices with nrows = length(contBC)
    pL2 = [pL2; (log(errL2(i,2:end))-log(errL2(i,1:end-1)))./(log(h(2:end))-log(h(1:end-1)))];
    pH1 = [pH1; (log(errH1(i,2:end))-log(errH1(i,1:end-1)))./(log(h(2:end))-log(h(1:end-1)))];
end


end
