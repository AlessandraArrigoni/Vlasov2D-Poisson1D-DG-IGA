function [adj_data] = get_adjoint_problem_data(pb_data, J_string, J_domain)
% Save the data for the adjoint problem starting from the original pde + BC
% and the considered functional (for now, we consider only functionals
% defined on the whole domain or on an full boundary)

neu_bdn_idx = pb_data.nmnn_sides; % only works when we have just ONE neumann boundary

adj_data.nmnn_sides = pb_data.nmnn_sides;
adj_data.drchlt_sides = pb_data.drchlt_sides;
adj_data.c_diff  = pb_data.c_diff; % works with the laplacian, to be checked if the diffusion is not 1
adj_data.grad_c_diff = pb_data.grad_c_diff;
adj_data.h = pb_data.h; % Dirichlet BC
        
switch J_domain
    case '0'
        switch J_string
            case 'u'
                adj_data.f = @(x, y) -1 + 0.*x.*y;
                adj_data.g = @(x, y, ind) 0.*x.*y;
                
            case 'u_squared'
                adj_data.f_type = '-2u_sol';
                adj_data.g = @(x, y, ind) 0.*x.*y;
                
            case 'grad_usquared'
                adj_data.f = @(x, y) -2*pb_data.f(x,y); 
                adj_data.g = @(x, y, ind) -2*pb_data.g(x,y);  
        end
                
    case string(neu_bdn_idx)
        switch J_string
            case 'u' % NOT WORKING IF NEUMANN BC != 1 (REDEFINE FUNCTION BELOW)
                adj_data.f = @(x, y) 0.*x.*y;
                adj_data.g = @neumannBC_functional;
                
            case 'u_squared'
                adj_data.f = @(x, y) zeros (size (x));
                adj_data.g_type  = '-2u_sol'; 
            
        end
        
    otherwise % no source and no neumann BC
        adj_data.f = @(x, y) 0.*x.*y;
        adj_data.g = @(x, y, ind) 0.*x.*y;    

end

end

function g = neumannBC_functional (x, y, ind)
  switch ind
    case 1
      g = -1 + 0.*x.*y;
    case 2
      g = 0.*x.*y;
    case 3
      g = 0.*x.*y;
    case 4
      g = 0.*x.*y;
    otherwise
      error ('g_nmnn: unknown reference number');
  end
end