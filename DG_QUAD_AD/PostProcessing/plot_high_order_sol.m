function [] = plot_high_order_sol (Dati, femregion, u_h, npoints1D)
% Plot the solution on each SQUARE element by evaluating the basis functions on
% a given set of points (not necessarily the quadrature nodes).
% We need this because we want to use high order polynomials and the matlab
% function patch connects the values in a linear way.

% INPUT : u_h is the set of coefficients with respect to the DG basis
%         npoints1D number of points in 1D where we want to evaluate the
%         solution (the first and last one are on the boundaries of the elem)

nln=femregion.nln;
ne=femregion.ne;

% points on the reference element [-1,1]^2
xref = linspace(-1,1,npoints1D);
points2D = []; % They have the same structure as node_2D from quadrature()
for i = 1:npoints1D
    for j = 1:npoints1D
        points2D = [points2D; xref(i), xref(j)];
    end
end


% scalar shape functions
[shape_basis]= basis_lagrange(Dati.fem);

% evaluation of shape functions
[dphiq, ~, ~]= evalshape(femregion, shape_basis, points2D, {xref,xref}, femregion.nln);
dphiq = reshape(dphiq,[length(points2D),femregion.nln]); % nqn2D x nln

% plot solution
figure()
x1=femregion.domain(1,1);
x2=femregion.domain(1,2);
y1=femregion.domain(2,1);
y2=femregion.domain(2,2);
M=max(u_h); % values for the z axis limits in the 3D plot
m=min(u_h);
if (abs(m-M) < 0.1)
    M=m+1;
end

for ie=1:ne % loop over elements
    
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    index_element=femregion.nedges*(ie-1).*ones(femregion.nedges,1) + [1:1:femregion.nedges]';

    coords_elem=femregion.coords_element(index_element, :);

    local_uh = u_h(index);
    
    sol_plot = zeros(npoints1D, npoints1D);
    
    for k = 1:length(points2D) % loop over evaluation points
        
        local_sol = dphiq(k,:)*local_uh;
        
        
        idx_i = mod(k-1, npoints1D) +1;
        idx_j = floor((k-1)/npoints1D) + 1;
        sol_plot(idx_i, idx_j) = local_sol;
        
    end
    
    % ADD PLOT OF THE SOLUTION (patch is good for linear functions)
    % translate the points on the current element and create the grid 
    [pphys_2D] = map_local_physical_points(coords_elem, points2D, femregion.BJ); %  Jacobian and physical coordinates at the quadrature points
    xx = pphys_2D(1 : npoints1D : end, 1);
    yy = pphys_2D(1 : npoints1D, 2);
    [XX, YY] = meshgrid(xx,yy);
    surf(XX,YY,sol_plot);
    shading('interp') % se metto 'faceted' plotta le linee e su ogni minielemento il colore è costante!
    hold on
    
    
end

axis([x1,x2,y1,y2,m,M]); colorbar;
%title(['Solution with fem = ' Dati.fem])
end