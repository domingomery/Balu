% ds      = Bcl_knn_old(X,d,Xt,options)  Training & Testing together
% options = Bcl_knn_old(X,d,options)     Training only
% ds      = Bcl_knn_old(Xt,options)      Testing only
%
% Toolbox: Balu
%    KNN (k-nearest neighbors) classifier. This implementation does not 
%    require VLFeat Toolbox. If you have it, please use Bcl_knn.
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.k is the number of neighbors (default=10)
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.string is a 8 character string that describes the performed
%       classification (e.g., 'knn,10  ' means knn with k=10).
%
%    Example: Training & Test together:
%       load datagauss                 % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)          % plot feature space
%       op.k = 10;
%       ds = Bcl_knn_old(X,d,Xt,op);   % knn with 10 neighbors
%       p = Bev_performance(ds,dt)     % performance on test data
%
%    Example: Training only
%       load datagauss                 % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)          % plot feature space
%       op.k = 10;
%       op = Bcl_knn_old(X,d,op);      % knn with 10 neighbors
%
%    Example: Testing only (after training only example):
%       ds = Bcl_knn_old(Xt,op);       % knn with 10 neighbors - testing
%       p = Bev_performance(ds,dt)     % performance on test data
%
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


function [ds,options] = Bcl_knn_old(varargin)

if exist('vl_kdtreequery','file')
    disp('Bcl_knn_old: You have installed VLFeat Toolbox. Please use Bcl_knn, it is faster.');
end

[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = sprintf('knn,%2d  ',options.k);
if train
    options.X      = X;
    options.d      = d;
    ds = options;
end
if test
    X = options.X;
    d = options.d;
    Nt = size(Xt,1);  
    N  = size(X,1);
    ds = zeros(Nt,1);
    k  = options.k;
    for q=1:Nt
       D = ones(N,1)*Xt(q,:)-X;
       [i,j] = sort(sum(D.*D,2));
       ds(q) = mode(d(j(1:k)));
    end
end
