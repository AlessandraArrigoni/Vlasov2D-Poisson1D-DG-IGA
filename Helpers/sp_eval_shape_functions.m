% SP_EVAL_SHAPE_FUNCTIONS: Compute the value or the derivatives of all the shape functions at a given set of points.
%
%   [val, F] = sp_eval_shape_functions ( space, geometry, pts, [options]);
%   [val, F] = sp_eval_shape_functions ( space, geometry, npts, [options]);
%
% INPUT:
%  
%     space:       object defining the discrete space (see sp_scalar)
%     geometry:    geometry structure (see geo_load)
%     pts:         cell array with coordinates of points along each parametric direction
%     npts:        number of points along each parametric direction
%     options:     cell array with the fields to plot
%                   accepted options are 'value' (default), 'gradient' and 'laplacian'
%
% OUTPUT:
%
%     val: cell array of matrices space.ndof x npts containing the values of all the
%          shape functions (or the derivatives) at all the points given as input
%     F:  grid points in the physical domain, that is, the mapped points
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2012, 2014, 2015, 2018 Rafael Vazquez
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

function [val, F] = sp_eval_shape_functions (space, geometry, npts, options)

  
  if (nargin < 4)
    options = {'value'};
  end
  if (~iscell (options))
    options = {options};
  end
  nopts = numel (options);

  ndim = numel (space.knots);
  
  endpoints = zeros (2, ndim);
  if (isfield (geometry, 'nurbs'))
    nurbs = geometry.nurbs;
    if (ndim == 1)
      nurbs.knots = {nurbs.knots};
    end
    for idim=1:ndim
      endpoints(:,idim) = nurbs.knots{idim}([nurbs.order(idim), end-nurbs.order(idim)+1]);
    end
    clear nurbs
  elseif (isfield (struct(space), 'knots'))
    degree = space.degree;
    for idim=1:ndim
      endpoints(:,idim) = space.knots{idim}([degree(idim)+1, end-degree(idim)]);
    end
  else
    endpoints(2,:) = 1;
  end
  
  if (iscell (npts))
    pts = npts;
    npts = cellfun (@numel, pts);
  elseif (isvector (npts))
    if (numel (npts) == 1)
      npts = npts * ones (1,ndim);
    end
    for idim = 1:ndim
      pts{idim} = linspace (endpoints(1,idim), endpoints(2,idim), npts(idim));
    end
  end

  for jj = 1:ndim
    pts{jj} = pts{jj}(:)';
    if (numel (pts{jj}) > 1)
      brk{jj} = [endpoints(1,jj), pts{jj}(1:end-1) + diff(pts{jj})/2, endpoints(2,jj)];
    else
      brk{jj} = endpoints(:,jj).';
    end
  end

  msh = msh_cartesian (brk, pts, [], geometry, 'boundary', false);
  sp  = space.constructor (msh);

  
  value = false; grad = false; laplacian = false; hessian = false;
  
  for iopt = 1:nopts
    switch (lower (options{iopt}))
      case 'value'
        val{iopt} = zeros(space.ndof, msh.nel); % assuming that the nel = npts
%         eunum{iopt} = {1:msh.nqn};
%         eusize{iopt} = npts;
        value = true;
        field{iopt} = 'shape_functions';

      case 'gradient'
        eu{iopt} = zeros (msh.rdim, msh.nqn, msh.nel);
%         eunum{iopt} = {1:msh.rdim, 1:msh.nqn};
%         eusize{iopt} = [msh.rdim, npts];
        grad = true;
        field{iopt} = 'shape_function_gradients';
%         
%       case 'laplacian'
%         eu{iopt} = zeros (msh.nqn, msh.nel);
%         eunum{iopt} = {1:msh.nqn};
%         eusize{iopt} = npts;
%         laplacian = true;
% 
%       case 'hessian'
%         eu{iopt} = zeros (msh.rdim, msh.rdim, msh.nqn, msh.nel);
%         eunum{iopt} = {1:msh.rdim, 1:msh.rdim, 1:msh.nqn};
%         eusize{iopt} = [msh.rdim, msh.rdim, npts];
%         hessian = true;
    end
  end

  F = zeros (msh.rdim, msh.nqn, msh.nel);

  msh_all = msh_precompute (msh);
  sp_all = sp_precompute (sp, msh_all, 'value', value, 'gradient', grad, ...
          'laplacian', laplacian, 'hessian', hessian);
      
% From here, it only works for 1D
    for iopt = 1:nopts
        switch (lower (options{iopt}))
            case 'value'
                sp_all.shape_functions = reshape (sp_all.shape_functions, sp_all.nsh_max, msh.nel);
            case 'gradient'
                sp_all.shape_function_gradients = reshape (sp_all.shape_function_gradients, sp_all.nsh_max, msh.nel);
        end
    end

  for iopt = 1:nopts
    for iel=1:msh.nel
        val{iopt}(sp_all.connectivity(:,iel), iel) = sp_all.(field{iopt})(:,iel);
    end
  end
  F = msh_all.geo_map;
  
%   for iel = 1:msh.nel_dir(1)
%     msh_col = msh_evaluate_col (msh, iel);
%     sp_col  = sp_evaluate_col (sp, msh_col, 'value', value, 'gradient', grad, ...
%           'laplacian', laplacian, 'hessian', hessian);
%     for iopt = 1:nopts
%      val{iopt}(sp_col.connectivity, iel) = sp_col.(field{iopt});    
%     end
%     F(:,:,msh_col.elem_list) = msh_col.geo_map;
% %     for iopt = 1:nopts
% %       eu{iopt}(eunum{iopt}{:},msh_col.elem_list) = eu_aux{iopt};
% %     end
%   end
  
  F = reshape (F, [msh.rdim, npts]);
%   for iopt = 1:nopts
%     eu{iopt} = reshape (eu{iopt}, [eusize{iopt}, 1]); % The extra 1 makes things work also in 1D
%   end

  if (nopts == 1)
    val = val{1};
  end
end
