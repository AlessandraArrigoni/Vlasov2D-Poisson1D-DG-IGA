% revert the entries of a vector

function y=swap_vector(x)

n=length(x);
y=[];

for k=1:1:n
    y(n-k+1)=x(k);
end