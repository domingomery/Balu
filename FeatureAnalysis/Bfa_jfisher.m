% J = Bfa_jfisher(X,d,p)
%
% Toolbox: Balu
%    Fisher objective function J.
%    X features matrix. X(i,j) is the feature j of sample i.
%    d vector that indicates the ideal classification of the samples
%    p a priori probability of each class
%
% See also Bfs_sfs.
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function J = Bfa_jfisher(X,d,p)

[n,M] = size(X);
dmin = min(d);

d = d-dmin+1;

N = max(d);

if not(exist('p','var'))
    p = ones(N,1)/N;
end

% Centroid of all samples
Xm = mean(X)';

L  = zeros(N,1);
Cw = zeros(M,M);
Cb = zeros(M,M);

for k=1:N
    ii   = find(d==k); % indices from class k
    L(k) = length(ii); % number of samples of class k
    Xk   = X(ii,:);    % samples of class k
    Xkm  = mean(Xk)';  % centroid of class k 
    Ck = cov(Xk);      % covariance of class k
    
    % within-class covariance
    Cw = Cw + p(k)*Ck; 
    
    % between-class covariance
    Cb = Cb + p(k)*(Xkm-Xm)*(Xkm-Xm)';
end

% Fisher discriminant
warning off
J = trace(Cw\Cb);
warning on
