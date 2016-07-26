% Sp = Bfa_sp100(X,d)
%
% Toolbox: Balu
%    Especificty at Sensibility = 100%.
%    X features matrix. X(i,j) is the feature j of sample i.
%    d vector that indicates the ideal classification of the samples
%
% See also Bfs_sfs, Bfa_fisher
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function Sp = Bfa_sp100(X,d)

[N,m] = size(X);

d = d-min(d);
if max(d)>1
    error('Bfa_sp100 works for only two classes');
end

d1 = find(d==1);

minz = min(X(d1,:));
maxz = max(X(d1,:));
z1 = zeros(size(X));
for p=1:m
    ii = find((X(:,p)>=minz(p))&(X(:,p)<=maxz(p)));
    z1(ii,p) = ones(length(ii),1);
end
if (m>1)
   drs = sum(z1,2);
   ii = find(drs==size(X,2));
   dr = zeros(size(d));
   dr(ii) = ones(length(ii),1);
else
   dr = z1;
end

TP = sum(dr.*d);
FP = sum(dr)-TP;

TN = sum(~dr.*~d);
Sp  = TN/(FP+TN);
