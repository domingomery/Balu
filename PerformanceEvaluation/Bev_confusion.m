% function [T,p] = Bev_confusion(d,ds,nn);
%
% Toolbox: Balu
%    Confusion Matrix and Performance of a classification
%
%    d is the ideal classification (vector Nx1 with N samples)
%    ds is the classified data (vector Nx1)
%    T is the confusion matrix (nxn) for n classes
%    T(i,j) indicates the number of samples i classified
%    as j.
%    p is the performance, diagonal sum divided by N
%    the classes should be labeled as 1, 2, ... n, but
%    if nn = [i j] is given it indicates that the lowest
%    class is min(nn) and the highest class is max(nn).
%    Thus, n = max(nn)-min(nn)+1 and the size of T is nxn.
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


function [T,p] = Bev_confusion(d,ds,nn)
if exist('nn','var');
    n1 = max(nn);
    n0 = min(nn);
else
    n1 = max([d;ds]);
    n0 = min([d;ds]);
end
n = n1-n0+1;
T = zeros(n,n);
for i=n0:n1
    for j=n0:n1
        kd = d==i;
        kds = ds==j;
        T(i-n0+1,j-n0+1) = sum(kd.*kds);
    end
end
p = trace(T)/sum(T(:));