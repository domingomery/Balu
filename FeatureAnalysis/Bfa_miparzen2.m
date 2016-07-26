% Mutual Information using Parzen windows for two variables
% NOTE:
%    The pdf's are estimated using Kernel Density Estimations programs 
%    kde.m and kde2d.m after Botev et al. (2010) implemented by Botev. 
%    These files are in Balu directory 'Feature Analysis' as Bfa_kde and
%    Bfs_kde2d. They can also be downloaded from www.mathwork.com
%    (c) Zdravko Botev. All rights reserved.

function I2 = Bfa_miparzen2(X,c,p)

[N,m] = size(X);
if m~=2
    error('This procedure is valid only for two features.');
end
if ~exist('c','var')
    c    = 0;
end

if ~exist('p','var')
    pvar    = 1;
else
    pvar    = 0;
end

%pxy = apdf2(X,T,1);
%px  = apdf1(X(:,1),T,0.02);
%py  = apdf1(X(:,2),T,0.02);

T = 256;

if c==0
    [bandwidth,pxy] = Bfa_kde2d(X,T);
    [bandwidth,px]  = Bfa_kde(X(:,1),T);
    [bandwidth,py]  = Bfa_kde(X(:,2),T);
else
    x = X(:,1);
    [bandwidth,px]  = Bfa_kde(x,T);
    d = X(:,2);
    xmin = min(x); xmax = max(x); dx = xmax-xmin;
    xmin = xmin-dx/10; xmax = xmax+dx/10;
    dmin = min(d); dmax = max(d);
    d = d-dmin+1; n = dmax-dmin+1;
    pxy = zeros(T,n);
    py  = zeros(n,1);
    for i=1:n
        ii = find(d==i);
        if pvar
            py(i) = p(i);
        else
            py(i) = length(ii)/N;
        end
        xi = x(ii);
        [bandwidth,pxi]  = Bfa_kde(xi,T,xmin,xmax);
        pxy(:,i) = pxi*py(i);
    end
end
pxy = pxy/sum2(pxy);
px  = px/sum(px);
py  = py/sum(py);

pxpy    = px*py';
t    = 1e-8;
I2   = sum(pxy(:).*log2((pxy(:)+t)./(pxpy(:)+t)));

