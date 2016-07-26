% function [ds,map] = Bct_medoidshift(X, sigma, k)
%
% Toolbox: Balu
%   
%    Medoidshift clustering
%    X matrix of samples.
%    sigma: standard deviation of the Gaussian Parzen window.
%    map is the map of the tree.
%    ds assigned class number.
%    k number of clusters.
%    map are the reduced coordinates.
%
%    Implementation based on:
%       Vedaldi,A; Stefano, S. (2008): Quick Shift and Kernel Methods for 
%       Mode Seeking, ECCV2008.
%
%    Example:
%       [X,d] = Bds_gaussgen([10 1;1 10;15 15],4*ones(3,3),100*ones(3,1));
%       figure(1)
%       Bio_plotfeatures(X,d);
%       ds = Bct_medoidshift(X,2,3);
%       figure(2)
%       Bio_plotfeatures(X,ds);
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [ds,map] = Bct_medoidshift(X, sigma, k)

G = X';
[d,N] = size(G) ;
oN = ones(N,1) ;
od = ones(d,1) ;
n = (G'.*G')*od ;
D = n*oN' + oN*n' - 2*(G'*G) ;
F = - exp(- .5 * D' / sigma^2) ;
Q = n * (oN'*F) - 2 * G' * (G*F) ;
[drop,map] = max(Q) ;
map = map';
ds = Bct_kmeans(map,k);
