% y = Bfa_gmean(X,d,op)
%
% Toolbox: Balu
%    op = 1: sqrt (Specificty * Sensibility)
%    op = 2: sqrt (Precision * Recall)
%    X features matrix. X(i,j) is the feature j of sample i.
%    d vector that indicates the ideal classification of the samples
%
% See also Bfs_sfs, Bfa_fisher
%
% (c) D.Mery, PUC-DCC, 2013
% http://dmery.ing.puc.cl

function y = Bfa_gmean(X,d,op)

[N,m] = size(X);

d = d-min(d);
if max(d)>1
    error('Bfa_gmean works for only two classes');
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
GT = sum

TN = sum(~dr.*~d);
Sp  = TN/(FP+TN);

switch op
    case 1
        y = sqrt()
    case 2
        PR = TP/(TP+FP);
        RE = TP/GT;
        


