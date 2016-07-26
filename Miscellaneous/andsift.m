% It separates sift descriptors into descriptors that belong to the region
% of interest R. R is a binary image. fi,di are the frames and descriptors
% of pixels = i (for i=0,1). (ff,dd) is the transposed output of vl_sift.

function [f1,d1,f0,d0,i1,i0] = andsift(ff,dd,R)

f = ff';
d = dd';

[N,M] = size(R);
ii    = round(f(2,:));
jj    = round(f(1,:));
kk = (ii<1)|(jj<1)|(ii>N)|(jj>M);
ii(kk) = [];
jj(kk) = [];

kk    = sub2ind([N M],ii,jj);
if ~isempty(kk)
    r1    = R(kk);
    i1    = find(r1==1);
    i0    = find(r1==0);
    f1    = f(:,i1);
    d1    = d(:,i1);
    f0    = f(:,i0);
    d0    = d(:,i0);
else
    f1 = [];
    d1 = [];
    f0 = f;
    d0 = d;
end
f1 = f1';
d1 = d1';
f0 = f0';
d0 = d0';