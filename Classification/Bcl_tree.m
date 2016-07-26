% ds      = Bcl_tree(X,d,Xt,options)  Training & Testing together
% options = Bcl_tree(X,d,options)     Training only
% ds      = Bcl_tree(Xt,options)      Testing only
%
% Toolbox: Balu
%    Classifier using a tree algorithm
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options is a strcuture containing a Balu classifier
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.n is the number of nodes
%       options.opt are the options of each classifier
%       options.im are the classes discriminated in each node
%       options.string is a 8 character string that describes the performed
%       classification (e.g., 'svm,4   ' means rbf-SVM).
%
%    Example: Training & Test together:
%       [X,d] = Bds_gaussgen([2 1;1 2;2 2;1 1],ones(4,2)/4,500*ones(4,1));
%       [Xt,dt] = Bds_gaussgen([2 1;1 2;2 2;1 1],ones(4,2)/4,500*ones(4,1));
%       Bio_plotfeatures(X,d)        % plot feature space
%       b.name = 'knn';
%       b.options.k=5;               % KNN with 5 neighbors
%       op = b;
%       ds = Bcl_tree(X,d,Xt,op);    % tree
%       p = Bev_performance(ds,dt)   % performance on test data
%
%    Example: Training only
%       [X,d] = Bds_gaussgen([2 1;1 2;2 2;1 1],ones(4,2)/4,500*ones(4,1));
%       [Xt,dt] = Bds_gaussgen([2 1;1 2;2 2;1 1],ones(4,2)/4,500*ones(4,1));
%       Bio_plotfeatures(X,d)        % plot feature space
%       b.name = 'knn';
%       b.options.k=5;               % KNN with 5 neighbors
%       op = b;
%       op = Bcl_tree(X,d,op);       % tree
%
%    Example: Testing only (after training only example):
%       ds = Bcl_tree(Xt,op);        % KNN with 5 neighbors
%       p = Bev_performance(ds,dt)   % performance on test data
%
%    Test the examples using a rbf-SVM classifier:
%       b.name = 'svm'; b.options.kernel = 4; op = b;
%       op = b;
%       op = Bcl_tree(X,d,op);       % tree with svm (training)
%       ds = Bcl_tree(Xt,op);        % tree with svm (test)
%       p = Bev_performance(ds,dt)   % performance on test data
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


function [ds,options] = Bcl_tree(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'tree    ';
if train
    d1 = min(d);
    d2 = max(d);
    n = d2-d1+1;

    XX = X;
    dd = d;

    im = zeros(n-1,1);
    t = 1;
    for k=n:-1:2
            perimax = -1;
            for i=1:n
                if sum(dd==i)>0
                    di = not(dd==i)+1;
                    [dsi,opi] = Bcl_structure(XX,di,XX,options);
                    peri = Bev_performance(dsi,di);
                    if peri>perimax
                        perimax = peri;
                        imax = i;
                        opx = opi;
                        opx.name = options.name;
                    end
                end
            end
            ii = find(dd==imax);
            XX(ii,:) = [];
            dd(ii)   = [];
            im(t) = imax;
            opt(t) = opx; %#ok<AGROW>
            t = t+1;
    end
    options.im = im;
    options.n   = n;
    options.opt = opt;
    ds = options;
end
if test
    XXt = Xt;
    im = options.im;
    n = options.n;
    ds = zeros(size(Xt,1),1);
    for i=1:n-1
        opi = options.opt(i);
        dsi = Bcl_structure(XXt,opi);
        ii = and(dsi==1,ds==0);
        ds(ii) = im(i);
    end
    ii = ds==0;
    h = (1:n)';
    h(im) = 0;
    jj = h>0;
    ds(ii) = h(jj);
end
