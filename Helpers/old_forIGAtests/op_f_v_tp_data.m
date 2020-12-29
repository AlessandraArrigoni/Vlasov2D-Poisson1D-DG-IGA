% OP_F_V_TP_DATA: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i), exploiting the tensor product structure.
%
%   rhs = op_f_v_tp_data (spv, msh, vett);
%
% INPUT:
%     
%   spv:   object representing the function space (see sp_scalar)
%   msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   vett:  vector of data containing the values of the rhs on the quadrature
%          nodes (we consider only the 1D case and they are ordered)
%
% OUTPUT:
%
%   rhs: assembled right-hand side
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

function rhs = op_f_v_tp_data (space, msh, vett)

  for idim = 1:msh.ndim
    size1 = size (space.sp_univ(idim).connectivity);
    if (size1(2) ~= msh.nel_dir(idim))
      error ('The discrete space is not associated to the mesh')
    end
  end
  
  rhs = zeros (space.ndof, 1);

  rhs = op_f_v(space, msh, vett);
%   for iel = 1:msh.nel_dir(1)
%     msh_col = msh_evaluate_col (msh, iel);
%     sp_col  = sp_evaluate_col (space, msh_col);
%     
% %     % Questo vale se vett è un vettore! 
% %     % Select the indices associated to the current element
% %     idx = (iel-1)*msh_col.nqn + 1 : iel*msh_col.nqn;
% %     rhs = rhs + op_f_v (sp_col, msh_col, vett(idx));
% 
%     % Questo vale se vett è una matrice con nrows=nqn in ogni elemento e
%     % ncols=nelementi 1D
%     rhs = rhs + op_f_v(sp_col, msh_col, vett(:,iel));
%   end

end
