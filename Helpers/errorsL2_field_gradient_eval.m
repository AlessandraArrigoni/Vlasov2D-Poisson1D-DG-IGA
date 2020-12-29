% errorsL2_field_gradient_eval: Evaluate the L2 error on a field living in 
% IGA space and on its gradient.
%
% NB : the input msh and space MUST HAVE ALREADY BEEN EVALUATED ! 
% INPUT:
%
%     space:    structure representing the space of discrete functions (see sp_scalar/sp_evaluate_col)
%     msh:      structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%     u:        vector of dof weights
%     uex:      function handle to evaluate the exact solution
%     graduex:  function handle to evaluate the gradient of the exact solution
%
% OUTPUT:
%
%     errL2:            error in L^2 norm on the field u
%     errL2_grad:       error in L^2 norm on the gradient \grad u
%     L2norm_exact:     L^2 norm of the exact solution uex
%     L2_norm_gradient: L^2 norm of the exact gradient graduex

function [errL2, errL2_grad, L2norm_exact, L2_norm_gradient] = errorsL2_field_gradient_eval (sp, msh, u, uex, graduex)

  grad_valu = sp_eval_msh (u, sp, msh, 'gradient');
  grad_valu = reshape (grad_valu, sp.ncomp, msh.rdim, msh.nqn, msh.nel);

  for idir = 1:msh.rdim
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end
  grad_valex  = reshape (feval (graduex, x{:}), sp.ncomp, msh.rdim, msh.nqn, msh.nel);

  w = msh.quad_weights .* msh.jacdet;

  [errL2, errl2_elem, L2norm_exact] = sp_l2_error (sp, msh, u, uex);
  % Error in H1 seminorm = L2 norm of the gradient
  errh1s_elem = sum (reshape (sum (sum ((grad_valu - grad_valex).^2, 1), 2), [msh.nqn, msh.nel]) .* w);
  L2norm_gradient_elem = sum (reshape (sum (sum (grad_valex.^2, 1), 2), [msh.nqn, msh.nel]) .* w);
 
  errL2_grad = sqrt (sum (errh1s_elem)); % Error in H1 seminorm = L2 norm of the gradient
  L2_norm_gradient = sqrt(sum(L2norm_gradient_elem));
  
   
end
