function [energy0, charge0, L2_norm0, ds, rhs_1, rhs_v2] = ...
        initialize_diagnostics(f0, E0, M_dg, msh, femregion, basis )
    
% PURPOSE: Compute the initial values for the L2 norm of f, the total charge and
% the discrete energy. We use the same number of quadrature nodes employed to
% compute the solutions (while for the errors we use twice as many nodes).
% Returns some vectors that are reused to compute the same quantities at the
% following timesteps
%
% INPUT: f0 -> coefficients of the initial density f with respect to DG basis
%               (vector with femregion.ndof elements)
%        E0 -> values of the initial electric field on the 1D quadrature nodes
%               (matrix n_quad_nodes_per_el X n_el)
%        M_dg -> mass matrix associated to the 2D DG basis
%        msh -> geoPDEs structure (run msh.precompute in advance)
%        femregion, basis
%
% OUTPUT: energy0 -> initial value of total energy (scalar)
%         charge0 -> initial value of total charge (scalar)
%         L2_norm0 -> initial value of L2 norm of density f (scalar)
%         ds -> 1D scaled quadrature weigths needed to compute the potential
%               energy at all the timesteps (matrix n_quad_nodes_per_el X n_el)
%         rhs_1 -> values of \int_{\omega} \phi dx dv where \phi is a DG basis
%               function (vector with femregion.ndof elements)
%         rhs_v2 -> values of \int_{\omega} \phi |v|^2 dx dv where \phi is a DG basis
%               function (vector with femregion.ndof elements)

% Compute POTENTIAL ENERGY: 0.5 * \int E_h .^2 dx
ds = msh.jacdet .* msh.quad_weights;
energy_pot0 = 0.5 * sum(sum((E0.^2).*ds));

% Compute KINETIC ENERGY: 0.5 * \int f_h * |v|^2 dx dv
rhs_v2 = source_rhs_given_f(femregion , @(x,y) abs(y).^2, basis);
energy_kin0 = 0.5 * (f0') * rhs_v2;
    
% Compute TOTAL ENERGY
energy0 = energy_kin0 + energy_pot0;

% Compute TOTAL CHARGE: \int f_h dx dv
rhs_1 = source_rhs_given_f(femregion, @(x,y) 1+0.*x.*y, basis);
charge0 = f0' * rhs_1;

% Compute L2 NORM OF THE DENSITY: sqrt(\int |f_h|^2 dx dv)
L2_norm0 = sqrt( f0'*M_dg*f0 );


end