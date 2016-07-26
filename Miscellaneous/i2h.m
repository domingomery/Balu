function y = i2h(x)
% conversion inhomogeneous to homogeneous coordinates
n = size(x,2);
y = [x;ones(1,n)];