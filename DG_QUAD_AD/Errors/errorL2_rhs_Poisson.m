function [error, error_rel, L2norm_exact] = errorL2_rhs_Poisson(Data, basis, u_h, femregion, region, L2norm_exact)
% Computes the L2 error on the rhs for Poisson equation (i.e. on the integral in v of
% the Vlasov density). Not used in the final tests of the project.
% 
% INPUT: Data -> struct with input data. Must contain a field "v_integr" with
%               a lambda function for the exact integral \int_{\Omega_v} f(x,v) dv 
%        u_h -> coefficients of the solution with respect to the DG basis
%        L2norm_exact -> L2 norm of the exact right hand side (computed if not
%               given as input
%
% OUTPUT: error -> L2 norm of the error 
%         error_rel -> relative L2 norm of the error
%         L2norm_exact -> L2 norm of the exact right hand side 

rhs_exact_f = @(x) Data.rho_0 - Data.v_integr(x);

[density_integralv, points_x] = compute_vertical_line_integral(basis, u_h, femregion, region);
rhs_num = Data.rho_0*ones(size(density_integralv)) - density_integralv;

% Plot
% result = reshape(rhs_num, [numel(points_x),1]);
% xplot = linspace(Data.domain(1,1), Data.domain(1,2),100);
% figure(); plot(points_x, result, '+-', xplot , rhs_exact_f(xplot));
% title('Line integral along y')

% Compute integral in x: rhs_num is a matrix [n_nodes1D, nel_1D]; points_x is a
% vector with n_nodes1D*nel_1D elements.

hx = region.x_points1D(2) - region.x_points1D(1); % each element in x has the same length by assumption
rhs_exact = reshape(rhs_exact_f(points_x), size(rhs_num));
if nargin < 6
    L2norm_exact = sum((basis.w_1D{1}*(rhs_exact.^2))*hx);
    L2norm_exact = sqrt(L2norm_exact);
end

error = sum((basis.w_1D{1}*((rhs_exact-rhs_num).^2))*hx);
error = sqrt(error);
error_rel = error/L2norm_exact;

end

