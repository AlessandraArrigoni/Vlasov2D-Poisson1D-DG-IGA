function [energy, charge, L2normVlasov, L2normE] = compute_diagnostics(f_cur, E_cur, M_dg, ds, rhs_1, rhs_v2)
% PURPOSE: Compute the current values for the L2 norm of f, the total charge,
% the discrete energy and the L2 norm of the electric field (needed for Landau
% damping test, only if required as output).
% We use the same number of quadrature nodes employed to
% compute the solutions (while for the errors we use twice as many nodes).
% 
% INPUT: f_cur -> coefficients of the current density f with respect to DG basis
%               (vector with femregion.ndof elements)
%        E_cur -> values of the current electric field on the 1D quadrature nodes
%               (matrix n_quad_nodes_per_el X n_el)
%        M_dg -> mass matrix associated to the 2D DG basis
%        ds -> 1D scaled quadrature weigths needed to compute the potential
%               energy at all the timesteps (matrix n_quad_nodes_per_el X n_el)
%        rhs_1 -> values of \int_{\omega} \phi dx dv where \phi is a DG basis
%               function (vector with femregion.ndof elements)
%        rhs_v2 -> values of \int_{\omega} \phi |v|^2 dx dv where \phi is a DG basis
%               function (vector with femregion.ndof elements)
%
% OUTPUT: all scalars

% Compute POTENTIAL ENERGY: 0.5 * \int E_h .^2 dx
energy_pot = 0.5 * sum(sum((E_cur.^2).*ds));

% Compute L2 NORM OF THE ELECTRIC FIELD: sqrt(\int E_h.^2 dx)
if nargout == 4
    L2normE = sqrt(2*energy_pot);
end

% Compute KINETIC ENERGY: 0.5 * \int f_h * |v|^2 dx dv
energy_kin = 0.5 * (f_cur') * rhs_v2;

% Compute TOTAL ENERGY
energy = energy_kin + energy_pot;

% Compute TOTAL CHARGE:  \int f_h dx dv
charge = f_cur' * rhs_1;

% Compute L2 NORM OF THE DENSITY: sqrt(\int |f_h|^2 dx dv)
L2normVlasov = sqrt( f_cur'*M_dg*f_cur );


end