% [X1,d1,X2,d2] = Bds_stratify(X,d,s)
%
% Toolbox: Balu
%
%    Data Stratification (without replacement)
%
%    input: (X,d) means features and ideal classification
%    Bds_stratify takes randomily a portion s (s between 0 and 1) of each class
%    from (X,d) to build (X1,d1). The samples not used in (X1,d1) are
%    stored in (X2,d2).
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [X1,d1,X2,d2,i1,i2] = Bds_stratify(X,d,s)


[i1,i2] = Bds_ixstratify(d,s);
X1 = X(i1,:); d1 = d(i1);
X2 = X(i2,:); d2 = d(i2);



% dmin = int8(min(d));
% dmax = int8(max(d));
% 
% 
% X1 = [];
% d1 = [];
% X2 = [];
% d2 = [];
% i1 = [];
% i2 = [];
% 
% for k=dmin:dmax
%     ik = find(d==k);
%     Xk = X(ik,:);
%     dk = d(ik);
%     nk = length(dk);
%     rk = rand(nk,1);
%     [i,j] = sort(rk);
%     sk = ceil(s*nk);
%     if (sk>0)
%         X1 = [X1; Xk(j(1:sk),:)];
%         d1 = [d1; dk(j(1:sk))];
%         i1 = [ii; ik(j(1:sk))];
%     end
%     if (sk<nk)
%         X2 = [X2; Xk(j(sk+1:nk),:)];
%         d2 = [d2; dk(j(sk+1:nk))];
%         i2 = [i2; ik(j(sk+1:nk))];
%     end
% end
