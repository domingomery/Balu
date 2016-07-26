% [i1,i2] = Bds_ixstratify(d,s)
%
% Toolbox: Balu
%
%    Data Stratification (without replacement)
%
%    input: (d,s) ideal classification and portion 
%    Bds_sixtratify takes randomily a portion s (s between 0 and 1) of each class
%    from d to build indices i1. The indices of not used samples are stored in i2.
%
% D.Mery, PUC-DCC, 2013
% http://dmery.ing.puc.cl

function [i1,i2] = Bds_ixstratify(d,s)

dmin = int8(min(d));
dmax = int8(max(d));

i1 = [];
i2 = [];

for k=dmin:dmax
    ik = find(d==k);
    dk = d(ik);
    nk = length(dk);
    rk = rand(nk,1);
    [i,j] = sort(rk);
    sk = ceil(s*nk);
    if (sk>0)
        i1 = [i1; ik(j(1:sk))];
    end
    if (sk<nk)
        i2 = [i2; ik(j(sk+1:nk))];
    end
end
