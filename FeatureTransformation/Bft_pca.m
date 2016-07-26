% [Y,lambda,A,Xs,mx] = Bft_pca(X,m)
%
% Toolbox: Balu
%    Principal component analysis
%    X is the matrix feature.
%    m number of selected components or the energy 0<m<=1 (in this case it
%    will be selected the first d principal components that fulfill the
%    condition sum(lambda(1:d))/sum(lambda) >= energy. See Example 2)
%    Y contains the m principal components of X
%    lambda contains the standard deviation of each principal component
%    A is the transformation matrix from X to Y
%    Xs is the reconstructed X from A and Y.
%
%  Example 1:
%     X = double(imread('cameraman.tif')); % 256x256 pixels
%     [Y,lambda,A,Xs] = Bft_pca(X,30);     % 30 principal components
%     figure(1);bar(lambda/lambda(1))
%     figure(2);imshow([X Xs],[])
%
%  Example 2:
%     X = double(imread('cameraman.tif')); % 256x256 pixels
%     [Y,lambda,A,Xs] = Bft_pca(X,0.90);   % 90% of energy
%     figure(1);bar(lambda/lambda(1))
%     figure(2);imshow([X Xs],[])
%
% (c) GRIMA, PUC-DCC, 2011
% http://grima.ing.puc.cl
%

function [Y,lambda,A,Xs,mx,B] = Bft_pca(X,d,options)

if nargin==2;
    options = d;
end


[N,L] = size(X);

if isfield(options,'m')
    m = options.m;
else
    m = options;
end

if and(m>0,m<1)
    energy = m;
    m = L;
else
    energy = 0;
end


mx = mean(X);
MX = ones(N,1)*mx;
X0 = X - MX;
Cx = cov(X);
[A,Cy] = eigsort(Cx);

lambda = diag(Cy);

if energy>0
    sumlambda = sum(lambda);
    energylam = zeros(L,1);
    for i=1:L
        energylam(i) = sum(lambda(1:i))/sumlambda;
    end
    ii = find(energylam>energy);
    m = ii(1);
end
        
B = A(:,1:m);
Y = X0*B; % Y = X0*A; Y = Y(:,1:m); % the first m components
Xs = Y*B' + MX;