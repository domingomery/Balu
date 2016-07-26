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

function [ds,Xc] = Bct_neighbor2D(X,th)

m = size(X,1);
nimg = size(X,2);

N = size(X,3);
ds = zeros(N,1);
c = 1;

ds(1) = c;

while(sum(ds>0)<N)
    i0 = find(ds==0);
    X0 = X(:,:,i0);
    n0 = length(i0);
    Xc = zeros(m,nimg,n0);
    for k=1:n0
        Xc(:,:,k) = mean(X0(:,:,k),2)*ones(1,nimg);
    end
    Xd = Xc-X0;
    M2 = sqrt(sum(Xd.*Xd,1));
    t = mean(M2,2);
    i2 = find(t<th);
    if ~isempty(i2)
        ds(i0(i2)) = c;
    else
        c = c+1;
        ds(i0(1))=c;
    end
end
Xc = zeros(m,nimg,c);

for i=1:c
    Xi = X(:,:,ds==i);
    ni = sum(ds==i);
    T   = zeros(nimg,ni);
    T(:) = sum(Xi,1);
    T = T';
    for p=1:m
        for q=1:nimg
            s = 0;
            nr = 0;
            for r=1:ni
                if T(r,q)>0
                    s = s + Xi(p,q,r);
                    nr = nr+1;
                end
            end
            if nr>0
               Xc(p,q,i) = s/nr;
            end
        end
    end
end
   