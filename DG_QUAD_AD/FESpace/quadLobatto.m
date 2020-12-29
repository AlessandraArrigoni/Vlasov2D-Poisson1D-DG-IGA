function [node_1D, w_1D, node_2D, w_2D] = quadLobatto(n, type_mesh)

if n == 3
    
    node_1D = [-1,0,1];
    w_1D_orig = [1/3, 4/3, 1/3];
    w_1D = w_1D_orig/2;
    
elseif n == 4
    
    node_1D = [-1,-sqrt(1/5), sqrt(1/5),1];
    w_1D_orig = [1/6, 5/6, 5/6, 1/6];
    w_1D = w_1D_orig/2;
    
elseif n == 5
    
    node_1D = [-1, -sqrt(3/7), 0, sqrt(3/7), 1];
    w_1D_orig = [1/10, 49/90, 32/45, 49/90, 1/10];
    w_1D = w_1D_orig/2;
    
elseif n == 6
    
    node_1D = [-1, -sqrt(1/3 + 2*sqrt(7)/21), -sqrt(1/3 - 2*sqrt(7)/21), sqrt(1/3 - 2*sqrt(7)/21),sqrt(1/3 + 2*sqrt(7)/21), 1];
    w_1D_orig = [1/15, (14-sqrt(7))/30,(14+sqrt(7))/30, (14+sqrt(7))/30,(14-sqrt(7))/30 ,1/15];
    w_1D = w_1D_orig/2;
    
end
       
        x = node_1D;
        w = w_1D_orig;
        node_2D=[];
        w_2D=[];
        for i=1:n
            for j=1:n
                node_2D=[node_2D; x(i), x(j)];
                w_2D = [w_2D, w(i).*w(j)];
            end
        end
end