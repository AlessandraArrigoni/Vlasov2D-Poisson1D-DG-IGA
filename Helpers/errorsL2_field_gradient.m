% PURPOSE: Compute the L2 norm of the errors and the L2 norms of the exact
% solution to a problem solved with IGA method.
%
% INPUT:
%
%    space:   object defining the space of discrete functions (see sp_scalar)
%    msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    graduex: function handle to evaluate the gradient of the exact solution
%
% OUTPUT:
%
%     errL2:            error in L^2 norm on the field u
%     errL2_grad:       error in L^2 norm on the gradient \grad u
%     L2norm_exact:     L^2 norm of the exact solution uex
%     L2_norm_gradient: L^2 norm of the exact gradient graduex


function [errL2, errL2_grad, L2norm_exact, L2_norm_gradient] = errorsL2_field_gradient (space, msh, u, uex, graduex)

  errL2 = 0; errL2_grad = 0;
  L2norm_exact = 0; L2_norm_gradient = 0;
  
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'value', true, 'gradient', true);
    
    [errL2_col, errL2_grad_col, L2norm_exact_col, L2_norm_gradient_col] = errorsL2_field_gradient_eval (sp_col, msh_col, u, uex, graduex);
    
    errL2_grad = errL2_grad + errL2_grad_col.^2;
    errL2 = errL2 + errL2_col.^2;
    L2norm_exact = L2norm_exact + L2norm_exact_col.^2;
    L2_norm_gradient = L2_norm_gradient + L2_norm_gradient_col.^2;
    
  end
  
  
  errL2 = sqrt (errL2);
  errL2_grad = sqrt (errL2_grad);
  L2norm_exact = sqrt(L2norm_exact);
  L2_norm_gradient = sqrt(L2_norm_gradient);

end
