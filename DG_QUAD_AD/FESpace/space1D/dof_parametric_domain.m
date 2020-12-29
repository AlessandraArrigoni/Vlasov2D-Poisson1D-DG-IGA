function [dofs] = dof_parametric_domain(femregion)
% Translates the ndof of the femregion 1D on the parametric domain (0,1) needed
% to evaluate the shape functions for IGA. I need to distinguish the various
% cases since the dof are not equispaced.
% I need this function since otherwise I get numerical approximation errors
% as I have [-pi,pi] as domain. 

dofs = zeros(femregion.ndof, 1);
coord_el = linspace(0,1,femregion.ne+1)'; % colonna
h = 1/femregion.ne;
nln = femregion.nln;

switch femregion.fem

    case 'P1'        
        for ie = 1:femregion.ne
            dofs((ie-1)*nln*ones(nln,1) + [1:nln]') = coord_el(ie:ie+1); 
        end
        
    case 'P2'
        for ie = 1:femregion.ne
            dofs((ie-1)*nln*ones(nln,1) + [1:nln]') = [coord_el(ie); coord_el(ie)+h/2; coord_el(ie+1)];
        end
        
    case 'P3'
        for ie = 1:femregion.ne
            dofs((ie-1)*nln*ones(nln,1) + [1:nln]') = [coord_el(ie); coord_el(ie)+h/4; coord_el(ie+1)-h/4; coord_el(ie+1)];
        end
end

end

