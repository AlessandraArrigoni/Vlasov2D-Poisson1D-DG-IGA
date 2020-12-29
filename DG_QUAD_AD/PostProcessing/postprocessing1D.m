function [] = postprocessing1D(femregion, data, uh, time, npoints, mytitle)
% PURPOSE: visualize the solution of a 1D linear transport problem on a
%          periodic domain, approximated by a DG polynomial space (P1, P2 or P3).
%          The exact solution u(x-trasp*time) given as lambda function in the
%          data struct is visualized too.
%
% INPUT : u_h -> set of coefficients with respect to the DG basis
%         npoints -> number of points in each 1D element where to evaluate the
%                    solution (the first and last one are the boundaries of the elem)

% Evaluate the shape functions on the given number of points
xx_ref = linspace(0,1,npoints); % basis functions are defined on (0,1)
shape1D = basis_lagrange1D(data.fem1D);
[phi1D, ~] = evalshape1D(shape1D, xx_ref, femregion.nln); 

% Plot numerical solution
xbeg = femregion.domain(1);
xend = femregion.domain(2);
M = max(uh); % values for the y axis limits
m = min(uh);
delta = M-m;
if (abs(m-M) < 0.1)
    M = m+1;
end


figure()
% Compute exact solution (be careful with the periodic domain!)
xx_plot = linspace(xbeg, xend, 100);
xx_periodic = xx_plot - data.trasp * time;
sol_exact = data.initial_f(xx_periodic,0); % The second entry is y coordinate

plot(xx_plot, sol_exact,'r--','LineWidth',1)
hold on
leg{1} = 'Exact';

% Numerical solution
for ie = 1:femregion.ne
    index = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]'; % indices of the 1D basis function on the current element
    coord_elem = femregion.coords_element((ie-1)*2*ones(2,1) + [1:2]'); % boundaries of the current element 
    
    % Select coefficients of the global solution linked to the current element
    uh_loc = uh(index);
    sol_loc = phi1D*uh_loc; % phi1D is a matrix npoints x nln
    
    % Set npoints on the current element
    xx_loc = coord_elem(1) + femregion.h*xx_ref; 
    
    plot(xx_loc, sol_loc,'b-','LineWidth',1)
    hold on   
    
end
leg{2} = 'Numerical';

legend(leg,'Location','NE','FontSize',14)
axis([xbeg,xend,m-delta/4,M+delta/4]); 
xlabel(femregion.direction);
if nargin > 5
    title(mytitle,'FontSize',18);
end

end