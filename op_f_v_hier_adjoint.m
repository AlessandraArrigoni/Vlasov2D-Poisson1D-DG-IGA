% OP_F_V_HIER_ADJOINT: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i) for hierarchical splines, 
%  exploiting the multilevel structure for the adjoint problem. It is the
%  same as the base one in the library but the input can be a spline
%  function (solution of the primal problem)
%
%   rhs = op_f_v_hier_adjoint (hspace, hmsh, f);
%   rhs = op_f_v_hier_adjoint (hspace, hmsh, u_sol, rhs_type);
%
% INPUT:
%     
%   hspace: object representing the hierarchical space of test functions (see hierarchical_space)
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   f:  function handle to compute the source function
%   u_sol: coefficients of the solution of the primal problem wrt hspace
%   rhs_type: string to define the rhs from the solution (should also work for Neumann boundary conditions)
%
% OUTPUT:
%
%   rhs: assembled right-hand side
% 
% The multilevel structure is exploited in such a way that only functions
%  of the same level of the active elements have to be computed. See also
%  op_gradu_gradv_hier for more details.
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function rhs = op_f_v_hier_adjoint (hspace, hmsh, varargin)
      
  rhs = zeros (hspace.ndof, 1);
  
  if (length(varargin) == 1) % function handle for the rhs --> assemble a vector 
      f = varargin{1};
      ndofs = 0;
  
      for ilev = 1:hmsh.nlevels
        ndofs = ndofs + hspace.ndof_per_level(ilev);
        if (hmsh.nel_per_level(ilev) > 0)
          x = cell (hmsh.rdim, 1);
          for idim = 1:hmsh.rdim
            x{idim} = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), hmsh.mesh_of_level(ilev).nqn, hmsh.nel_per_level(ilev));
          end
          sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
          sp_lev = change_connectivity_localized_Csub (sp_lev, hspace, ilev);
          b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, f(x{:}));

          dofs = 1:ndofs;
          rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}' * b_lev;
        end
      end
  elseif (length(varargin) == 2) % vector of coefficients and string to choose the type of rhs
      
      u_sol = varargin{1};
      f_type = varargin{2};
      
      if (f_type == '-2u_sol')
          mass_mat = op_u_v_hier(hspace, hspace, hmsh);
          rhs = -2 * mass_mat * u_sol;
      end
      
      % TODO : add cases where only a part of the domain or the boundary is
      % considered.
  else
      disp('ERROR: wrong number of input arguments');
  end

end
