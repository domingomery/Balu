% function [ds,Z] = Bct_meanshift(X, sigma)
%
% Toolbox: Balu
%   
%    Meanshift clustering
%    X matrix of samples
%    sigma: standard deviation of the Gaussian Parzen window
%    ds assigned class number.
%    Z are the reduced coordinates.
%
%    Implementation based on:
%       Vedaldi,A; Stefano, S. (2008): Quick Shift and Kernel Methods for 
%       Mode Seeking, ECCV2008.
%
%    Example:
%       [X,d] = Bds_gaussgen([10 1;1 10;15 15],4*ones(3,3),100*ones(3,1));
%       figure(1)
%       Bio_plotfeatures(X,d);
%       ds = Bct_meanshift(X,2);
%       figure(2)
%       Bio_plotfeatures(X,ds);
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [ds,Z] = Bct_meanshift(X, sigma)
G = X';
[d,N] = size(G) ;
oN = ones(N,1) ;
od = ones(d,1) ;
n = (G'.*G')*od ;
Z = G ;
T = 100 ;
for t=1:T
    m = (Z'.*Z')*od ;
    D = m*oN' + oN*n' - 2*(Z'*G) ;
    F = - exp(- .5 * D' / sigma^2) ;
    Y = F ./ (oN * (oN'*F)) ;
    Z = G*Y ;
end
Z = Z';
ds = zeros(N,1);
ds(1)=1;
ms = Z(1,:);
ns = 1;
r = 1;
s = sigma*2;
for i=2:N
    mi = Z(i,:);
    D = ms-ones(r,1)*mi;
    D2 = D.*D;
    Dk = sum(D2,2);
    [ii,jj] = min(Dk');
    if ii<s
        ds(i) = jj(1);
        ms(jj,:) = ms(jj,:)*ns(jj)+mi;
        ns(jj)=ns(jj)+1;
        ms(jj,:) = ms(jj,:)/ns(jj);
    else
        r = r+1;
        ms = [ms;mi];
        ns = [ns;1];
        ds(i)=r;
    end
end


