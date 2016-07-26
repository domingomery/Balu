function [Y,j] = Bft_vq(X,Xcen,kd)

if exist('kd','var')
    j       = vl_kdtreequery(kd,Xcen',X','NumNeighbors',1)';
else
    [jj, j] = min(vl_alldist(Xcen', X'), [], 1) ;
end

Y    = Xcen(j,:);
