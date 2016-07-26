% function [ds,Xc] = Bct_neighbor(X,th)
%
% Toolbox: Balu
%   
%    Neigbor clustering: iterative method, a sample will be added to a 
%    cluster if its distance to the mass center of the cluster is less than 
%    th, else it will be create a new cluster.
%
%    X matrix of samples
%    th minimal distance
%    ds assigned class number
%    Xc mass center of each cluster
%
%    Example:
%     [X,d] = Bds_gaussgen([10 1;1 10],1*ones(2,2),100*ones(2,1));
%     figure(1)
%     Bio_plotfeatures(X,d);
%     ds = Bct_neighbor(X,4);
%     figure(2)
%     Bio_plotfeatures(X,ds);
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [ds,Xc] = Bct_neighbor(X,th)

m = size(X,2);

N = size(X,1);
ds = zeros(N,1);
c = 1;

ds(1) = c;


while(sum(ds>0)<N)
    i0 = find(ds==0);
    X0 = X(i0,:);
    n0 = length(i0);
    Xc = ones(n0,1)*mean(X(ds==c,:),1);
    Xd = Xc-X0;
    M2 = sqrt(sum(Xd.*Xd,2));
    i2 = find(M2<th);
    if ~isempty(i2)
        ds(i0(i2)) = c;
    else
        c = c+1;
        ds(i0(1))=c;
    end
end
Xc = zeros(c,m);
for i=1:c
    Xc(i,:) = mean(X(ds==i,:),1);
end