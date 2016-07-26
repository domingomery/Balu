% [per,R] = Bfa_corrsearch(x,y,method,v,show)
%
% Toolbox: Balu
%    Error estimation of linear or quadratic model
%    that minimizes the norm between measured and modeled
%    output. The error is estimated using cross-validation.
%
%    x: measured input (nxm : n samples and m variables)
%    y: measured output (n measurements)
%    v: number of groups in cross-validation
%    method: 1 linear, 2 quadratic
%
% D.Mery, PUC-DCC, Aug. 2008
% http://dmery.ing.puc.cl
%
function [per,R] = Bfa_corrsearch(x,y,method,v,show)

[n,m] = size(x);
switch method
    case 1 % linear model
        X = [ones(n,1) x];
    case 2 % quadratic model
        p = [nchoosek(1:m,2);[(1:m)' (1:m)']];
        X = [ones(n,1) x];
        for i=1:size(p,1)
            X = [X x(:,p(i,1)).*x(:,p(i,2))];
        end
end


if not(exist('v','var'))
    v = 10;
end

if (v==1)
    error('cross validation does not work with only one group.');
end


N = size(X,1);

rn = rand(N,1);
[i,j] = sort(rn);

Xr = X(j,:);
dr = y(j);

r = fix(N/v);
R = zeros(v,2);
ini = 1;
for i=1:v-1
    R(i,:) = [ini ini+r-1];
    ini = ini + r;
end
R(v,:) = [ini N];

per = [];

for i=1:v
    XXt = Xr(R(i,1):R(i,2),:);
    ddt = dr(R(i,1):R(i,2),:);
    XX = [];
    dd = [];
    for j=1:v
        if (j~=i)
            XX = [XX;Xr(R(j,1):R(j,2),:)];
            dd = [dd;dr(R(j,1):R(j,2),:)];
        end
    end

    a = inv(XX'*XX)*XX'*dd;
    dds = XXt*a;
    pp = mean(abs(dds-ddt));
    per = [per; pp;];
end

% confidence intervals for c = 95%
per = per;
c = 0.95;
p     = mean(per);
mu    = p;
sigma = sqrt(p*(1-p)/v);
t = (1-c)/2;
if v==size(X,1)
    z = norminv(1-t);
    p1 = mu - z*sigma;
    p2 = mu + z*sigma;
else
    z = tinv(1-t,v-1);
    vv = sigma^2;
    p1 = mu - z*vv/sqrt(v);
    p2 = mu + z*vv/sqrt(v);
end
if show
    disp(sprintf('%7.4f < ths = %7.4f < %7.4f with %3.0f%% confidence',p1*100,p*100,p2*100,c*100));
end
a = inv(X'*X)*X'*y;
ys = X*a;

X = [y ones(n,1)];
a = inv(X'*X)*X'*ys;

R = corr(y,ys);
%subplot(2,2,kl)
if (show)
    plot(y,ys,'.')
    xlabel('measured')
    ylabel('modeled')
    ax = axis;
    XX = [ax(1:2)' ones(2,1)];
    YY = XX*a;
    hold on
    plot(XX(:,1),YY,'r')
    text(0.5,5,sprintf('R=%5.4f',R))
    text(0.5,4.5,sprintf('%5.4f < e = %5.4f < %5.4f',p1,p,p2));
    %title(tl)
end