%--------------------------------------------------------------------
% PURPOSE:
%
% This routine makes the post-processing by plotting the computational 
% mesh as well as the exact and the discrete solutions.
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [solutions]= postprocessing(femregion,Dati,u_h)

x=femregion.dof(:,1);
y=femregion.dof(:,2);

if isfield(Dati, 'exact_sol')
    u_ex = Dati.exact_sol(x,y,Dati.time);
else
    u_ex = Dati.initial_f(x,y);
end

% % plot solution
% figure
% stem3(femregion.dof(:,1),femregion.dof(:,2),u_ex,'-r*');
% hold on
% stem3(femregion.dof(:,1),femregion.dof(:,2),u_h,'-bo');
% legend('Exact Solution', 'Computed solution')
% hold off

% plot solution
figure()
x1=femregion.domain(1,1);
x2=femregion.domain(1,2);
y1=femregion.domain(2,1);
y2=femregion.domain(2,2);
M=max(u_h);
m=min(u_h);
if (abs(m-M) < 0.1)
    M=m+1;
end


if Dati.fem == 'Q1'
    k = 1;
    for ie = 1:femregion.ne
        patch(femregion.dof(k:k+3,1), femregion.dof(k:k+3,2), full(u_h(k:k+3)), full(u_h(k:k+3)));
        view(3)
        hold on
        k = k+4;
    end    

elseif Dati.fem == 'Q2'
    k = 1;
    for ie = 1 : femregion.ne
        % The points are the vertices of each square in counterclockwise order
        % starting from top left
        patch(femregion.dof([k,k+2,k+8,k+6],1), femregion.dof([k,k+2,k+8,k+6],2), full(u_h([k,k+2,k+8,k+6])), full(u_h([k,k+2,k+8,k+6])) ); %patch(x, y, values, colors)
        view(3)
        hold on;
        k=k+9;
    end   
    
elseif Dati.fem == 'Q3'
    k = 1;
    for ie = 1 : femregion.ne
        % The points are the vertices of each square in counterclockwise order
        % starting from top left
        patch(femregion.dof([k,k+3,k+12,k+15],1), femregion.dof([k,k+3,k+12,k+15],2), full(u_h([k,k+3,k+12,k+15])), full(u_h([k,k+3,k+12,k+15])) ); %patch(x, y, values, colors)
        view(3)
        hold on;
        k=k+16;
    end   
    
elseif Dati.fem == 'P1'    
    k = 1;
    for ie = 1 : femregion.ne
        trisurf([1 2 3],femregion.dof(k:k+2,1),femregion.dof(k:k+2,2),full(u_h(k:k+2)))
        hold on;
        k=k+3;
    end
    
    
elseif Dati.fem == 'P2'
    k = 1;
    for ie = 1 : femregion.ne
        trisurf([1 2 3],femregion.dof([k,k+2,k+4],1),femregion.dof([k,k+2,k+4],2),full(u_h([k,k+2,k+4])))
        hold on;
        k=k+6;
    end   
    
elseif Dati.fem == 'P3'
    k = 1;
    for ie = 1 : femregion.ne
        trisurf([1 2 3],femregion.dof([k,k+3,k+6],1),femregion.dof([k,k+3,k+6],2),full(u_h([k,k+3,k+6])))
        hold on;
        k=k+10;
    end       
end

view(2)
title(['u_h(x,y) with fem = ' Dati.fem ', nref = ' num2str(femregion.nref)]); xlabel('x-axis'); ylabel('y-axis');
axis([x1,x2,y1,y2,m,M]); colorbar;

solutions=struct('u_h',u_h,...
    'u_ex',u_ex);
end
