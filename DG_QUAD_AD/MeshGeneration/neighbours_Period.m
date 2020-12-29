%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create TOPOLOGICAL MATRIX of the neighbours: for each element we specify its
% neighbours along each of the sides.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines neighboring elements and number of neighboring edges of
% a list of elements given by connectivity list t
%
% neigh(i,j) = jth neighbor of element i (the order: left, bottom, right, top)
% neighedges(i,j) = number of edge of jth neighbor of element i
%
% 1 \leq i \leq number of elements
% 1 \leq j \leq nedges
%
% value: -1 if boundary edge. We assume that the mesh is PERIODIC in x, so we
% can have boundary edges only at the top or bottom of the elements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [neighbor]= neighbours_Period(femregion)

fprintf('============================================\n')
fprintf('LOG:: Getting neighbors\n')
fprintf('============================================\n')

switch femregion.type_mesh

    case{'TS', 'TU'}

        NEL=femregion.ne;
        connectivity=femregion.connectivity;

        neigh=-ones(NEL,3);
        neighedges=-ones(NEL,3);

        for i=1:(NEL-1)

            v1=connectivity(1,i);
            v2=connectivity(2,i);
            v3=connectivity(3,i);

            for j=(i+1):NEL

                vn1=connectivity(1,j);
                vn2=connectivity(2,j);
                vn3=connectivity(3,j);

                flag=0;
                if (v1==vn2 & v2==vn1)
                    neigh(i,1)=j;
                    neigh(j,1)=i;
                    neighedges(i,1)=1;
                    neighedges(j,1)=1;
                elseif (v1==vn1 & v2==vn3)
                    neigh(i,1)=j;
                    neigh(j,3)=i;
                    neighedges(i,1)=3;
                    neighedges(j,3)=1;
                elseif (v1==vn3 & v2==vn2)
                    neigh(i,1)=j;
                    neigh(j,2)=i;
                    neighedges(i,1)=2;
                    neighedges(j,2)=1;
                elseif (v2==vn2 & v3==vn1)
                    neigh(i,2)=j;
                    neigh(j,1)=i;
                    neighedges(i,2)=1;
                    neighedges(j,1)=2;
                elseif (v2==vn1 & v3==vn3)
                    neigh(i,2)=j;
                    neigh(j,3)=i;
                    neighedges(i,2)=3;
                    neighedges(j,3)=2;
                elseif (v2==vn3 & v3==vn2)
                    neigh(i,2)=j;
                    neigh(j,2)=i;
                    neighedges(i,2)=2;
                    neighedges(j,2)=2;
                elseif (v3==vn2 & v1==vn1)
                    neigh(i,3)=j;
                    neigh(j,1)=i;
                    neighedges(i,3)=1;
                    neighedges(j,1)=3;
                elseif (v3==vn1 & v1==vn3)
                    neigh(i,3)=j;
                    neigh(j,3)=i;
                    neighedges(i,3)=3;
                    neighedges(j,3)=3;
                elseif (v3==vn3 & v1==vn2)
                    neigh(i,3)=j;
                    neigh(j,2)=i;
                    neighedges(i,3)=2;
                    neighedges(j,2)=3;
                end

            end
        end

    case{'CART'}
        % Perform special checks for the elements on the boundary to enforce
        % periodic BC along the x direction; we must check both sides, not
        % only the left one since the loop over j starts from i+1 (we don't
        % consider all the couples each time, so we may miss something if the
        % element on the left boundary comes "after" the corresponding one on
        % the right border)

        ne=femregion.ne; 
        nedge = femregion.nedges; % =4
        
        % Initialize neighbour matrix
        neigh=-ones(ne,nedge);
        neighedges = -ones(ne,nedge);

        NEL=femregion.ne;
        connectivity=femregion.connectivity;
        xmin = femregion.domain(1,1);
        xmax = femregion.domain(1,2);

        for i=1:(NEL-1)

            v1=connectivity(1,i);
            v2=connectivity(2,i);
            v3=connectivity(3,i);
            v4=connectivity(4,i);
            
            % Get coordinates of the vertices
            c1 = femregion.coords_element(nedge*(i-1)+1,:);
            c2 = femregion.coords_element(nedge*(i-1)+2,:);
            c3 = femregion.coords_element(nedge*(i-1)+3,:);
            c4 = femregion.coords_element(nedge*(i-1)+4,:);

            for j=(i+1):NEL

                vn1=connectivity(1,j);
                vn2=connectivity(2,j);
                vn3=connectivity(3,j);
                vn4=connectivity(4,j);
                
                % Get coordinates of the vertices
                cn1 = femregion.coords_element(nedge*(j-1)+1,:);
                cn2 = femregion.coords_element(nedge*(j-1)+2,:);
                cn3 = femregion.coords_element(nedge*(j-1)+3,:);
                cn4 = femregion.coords_element(nedge*(j-1)+4,:);

                % Neighbours along the horizontal direction (periodic BC)
                if ((v1==vn4 && v2==vn3) || ...
                     (c1(1)==xmin && cn4(1)==xmax && c1(2)==cn4(2) && ...
                       c2(1)==xmin && cn3(1)==xmax && c2(2)==cn3(2)) || ...
                     (c1(1)==xmax && cn4(1)==xmin && c1(2)==cn4(2) && ...
                       c2(1)==xmax && cn3(1)==xmin && c2(2)==cn3(2)))
                    neigh(i,1)=j;
                    neigh(j,3)=i;
                    neighedges(i,1)=3;
                    neighedges(j,3)=1;

                elseif ((v3==vn2 && v4==vn1) || ...
                        (c3(1)==xmax && cn2(1)==xmax && c3(2)==cn2(2) && ...
                            c4(1)==xmin && cn1(1)==xmax && c4(2)==cn1(2)) || ...
                        (c3(1)==xmax && cn2(1)==xmin && c3(2)==cn2(2) && ...
                            c4(1)==xmax && cn1(1)==xmin && c4(2)==cn1(2)))
                    neigh(i,3)=j;
                    neigh(j,1)=i;
                    neighedges(i,3)=1;
                    neighedges(j,1)=3;

                % Neighbours along the vertical direction
                elseif (v2==vn1 & v3==vn4)
                    neigh(i,2)=j;
                    neigh(j,4)=i;
                    neighedges(i,2)=4;
                    neighedges(j,4)=2;



                elseif (v1==vn2 & v4==vn3)
                    neigh(i,4)=j;
                    neigh(j,2)=i;
                    neighedges(i,4)=2;
                    neighedges(j,2)=4;
                end
            end
        end

    otherwise

        error('The mesh could only be of triangular elements or quadrilateral elements');
end

neighbor=struct( 'neigh',neigh,...
    'neighedges',neighedges,...
    'nedges',femregion.nedges);

end

