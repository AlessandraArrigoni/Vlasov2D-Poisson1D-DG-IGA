%--------------------------------------------------------------------
% PURPOSE:
%
%  This routine plots a PDE triangular/quadrilateral mesh.
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------
function plot_mesh(vertices, connectivity, domain)  

% fill the domain
fill([domain(1,1) domain(1,2) domain(1,2) domain(1,1) domain(1,1)], ...
     [domain(2,1) domain(2,1) domain(2,2) domain(2,2) domain(2,1)], ...
     'w','LineWidth',2.0)    


[r,c] = size(connectivity);

hold on;

for ie=1:c
    
    index_el=(r-1)*(ie-1).*ones(1,r-1) + [1:1:(r-1)];
    
    if r==4
        fill([vertices(connectivity(1,ie),1) vertices(connectivity(2,ie),1) ...
                vertices(connectivity(3,ie),1) vertices(connectivity(1,ie),1)], ...
            [vertices(connectivity(1,ie),2) vertices(connectivity(2,ie),2) ...
                vertices(connectivity(3,ie),2) vertices(connectivity(1,ie),2)],'w','LineWidth',2.0)
        
    else
        
        fill([vertices(connectivity(1,ie),1) vertices(connectivity(2,ie),1) ...
                vertices(connectivity(3,ie),1) vertices(connectivity(4,ie),1) ...
                vertices(connectivity(1,ie),1)], ...
            [vertices(connectivity(1,ie),2) vertices(connectivity(2,ie),2) ...
                vertices(connectivity(3,ie),2) vertices(connectivity(4,ie),2) ...
                vertices(connectivity(1,ie),2)],'w','LineWidth',2.0)
    end
end

axis square
axis tight
set(gca,'XTick',[])
set(gca,'XTickLabel','')
set(gca,'YTick',[])
set(gca,'YTickLabel','')
