function y = h2i(x)
% conversion homogeneous to inhomogeneous coordinates
n = size(x,1);
y = x(1:n-1,:)./(ones(n-1,1)*x(n,:));