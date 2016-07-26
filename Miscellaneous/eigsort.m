function [A,D] = eigsort(X)
L = size(X,2);
[Ao,Do] = eig(X);
A = Ao(:,L:-1:1);
D = Do(:,L:-1:1);
D = D(L:-1:1,:);
