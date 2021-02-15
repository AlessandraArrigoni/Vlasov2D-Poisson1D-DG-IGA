function [value] = compute_functionalJ(Jstring, Jdomain, u, hmsh, hspace)

% COMPUTE_FUNCTIONALJ : Compute the value of the functional J on the current domain,
% starting from the overkilled solution u in hspace. 
% We use the function "sp_h1_error" with the fake exact solution u_ex = 0. 
% !!! For now only functionals defined on the WHOLE domain or on one FULL
% boundary are considered !!!
%
% INPUT: 
%           Jstring: 'u_squared','grad_usquared'; for now we don't consider
%                   the case 'u' otherwise we can't use sp_h1_error and we
%                   should rewrite the functions.
%           Jdomain:  0 = interior, 1,2,3,4 = left, right, bottom, top
%                   boundary; for now we cannot compute gradients on the boundary
%                   (we don't know the adjoint problem)
%           u:       solution coefficients in the overkilled THB-space
%           hmsh/hspace: hierarchical mesh and space for the overkilled
%                    solution u. hmsh stores the mapping for the current
%                    geometry.
%
% OUTPUT:
%           value: functional value
%
%


switch Jdomain
    case '0' % interior
        switch Jstring
            case 'u'
                disp('ERROR: not implemented functional: \int_{\Omega}(u)'); 
            case 'u_squared'
                [value, ~] = sp_l2_error (hspace, hmsh, u, @(x,y) 0.*x.*y);

            case 'grad_usquared' % seminorm of the H1 error
                graduex = @(x, y) cat (1, ...
                       reshape (0.*x.*y, [1, size(x)]), ...
                       reshape (0.*x.*y, [1, size(x)]));
                [~, ~, value] = sp_h1_error(hspace, hmsh, u, @(x,y) 0.*x.*y, graduex);
        end        
        
    otherwise
        side = str2double(Jdomain);
        switch Jstring
            case 'u'
                disp('ERROR: not implemented functional: \int_{\Gamma}(u)');
            case 'u_squared'
                u_bdn = u(hspace.boundary(side).dofs);
                [value, ~] = sp_l2_error (hspace.boundary(side), hmsh.boundary(side), u_bdn, @(x,y) 0.*x.*y); % we pass uex as a function of 2 variables since the hmsh.boundary has 2 components
                
            case 'grad_usquared'
                disp('ERROR: not implemented functional: \int_{\Gamma}(grad_u^2)'); 
        end       
        
        
end



end

