% ds      = Bcl_knn(X,d,Xt,options)  Training & Testing together
% options = Bcl_knn(X,d,options)     Training only
% ds      = Bcl_knn(Xt,options)      Testing only
%
% Toolbox: Balu
%    KNN (k-nearest neighbors) classifier using randomized kd-tree
%    forest from FLANN. This implementation requires VLFeat Toolbox.
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
%       options.kdtree contains information about the randomized kdtree
%       (from function vl_kdtreebuilf of VLFeat Toolbox).
%       options.string is a 8 character string that describes the performed
%       classification (e.g., 'knn,10  ' means knn with k=10).
%
%    Example: Training & Test together:
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.k = 10;
%       ds = Bcl_knn(X,d,Xt,op);   % knn with 10 neighbors
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.k = 10;
%       op = Bcl_knn(X,d,op);      % knn with 10 neighbors
%
%    Example: Testing only (after training only example):
%       ds = Bcl_knn(Xt,op);       % knn with 10 neighbors - testing
%       p = Bev_performance(ds,dt) % performance on test data
%
%
% D.Mery, C. Mena PUC-DCC, 2010-2013
% http://dmery.ing.puc.cl


function [ds,options] = Bcl_knn(varargin)

if ~exist('vl_kdtreequery','file')
    error('Bcl_knn: This function requires the VLFeat Toolbox.');
end

[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = sprintf('knn,%2d  ',options.k);
if train
    options.kdtree = vl_kdtreebuild(X');
    options.X      = X;
    options.d      = d;
    ds = options;
end
if test
    [i,dist] = vl_kdtreequery(options.kdtree,options.X',Xt','NumNeighbors',options.k);
    if (options.k > 1)
        ds       = mode(options.d(i))';
    else % for 1-KNN (modification by Carlos Mena)
        ds       = options.d(i);
    end
    
    if isfield(options,'output')
        ns = length(ds);
        sc = zeros(ns,1);
        for q=1:ns
            j = find(options.d(i(:,q))==ds(q));
            sc(q) = min(dist(j,q));
        end
        ds = Bcl_outscore(ds,sc,options);
    end
end
