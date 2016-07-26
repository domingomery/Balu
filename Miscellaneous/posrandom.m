% It changes the order of the columns or rows of x randomly. Variable dim 
% selects the dimension along which to change.

function [x_new,j] = posrandom(x,dim)

nx    = size(x,dim);
r     = rand(nx,1);
[i,j] = sort(r);
if dim==1
    x_new = x(j,:);
else
    x_new = x(:,j);
end
