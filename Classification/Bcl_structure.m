% ds      = Bcl_structure(X,d,Xt,options)  Training & Testing together
% options = Bcl_structure(X,d,options)     Training only
% ds      = Bcl_structure(Xt,options)      Testing only
%
% Toolbox: Balu
%    Classification using Balu classifier(s) defined in structure b.
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options is a Balu classifier structure b with
%          b.name      = Balu classifier's name
%          b.options   = options of the classifier
%
%       b can define one or more classifiers (see example).
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data (one column per classifier)
%
%    Example: Training & Test together:
%       load datagauss                                                        % simulated data (2 classes, 2 features)
%       b(1).name = 'knn';   b(1).options.k = 5;                              % KNN with 5 neighbors
%       b(2).name = 'knn';   b(2).options.k = 7;                              % KNN with 7 neighbors
%       b(3).name = 'knn';   b(3).options.k = 9;                              % KNN with 9 neighbors
%       b(4).name = 'lda';   b(4).options.p = [];                             % LDA
%       b(5).name = 'qda';   b(5).options.p = [];                             % QDA
%       b(6).name = 'nnglm'; b(6).options.method = 3; b(6).options.iter = 10; % Nueral network
%       b(7).name = 'svm';   b(7).options.kernel = 4;                         % rbf-SVM
%       b(8).name = 'maha';  b(8).options = [];                               % Euclidean distance
%       b(9).name = 'dmin';  b(9).options = [];                               % Mahalanobis distance
%       op = b;
%       ds = Bcl_structure(X,d,Xt,op);                                        % ds has 9 columns
%       p = Bev_performance(ds,dt)                                            % p has 9 performances
%
%
%    Example: Training only
%       load datagauss                                                        % simulated data (2 classes, 2 features)
%       b(1).name = 'knn';   b(1).options.k = 5;                              % KNN with 5 neighbors
%       b(2).name = 'knn';   b(2).options.k = 7;                              % KNN with 7 neighbors
%       b(3).name = 'knn';   b(3).options.k = 9;                              % KNN with 9 neighbors
%       b(4).name = 'lda';   b(4).options.p = [];                             % LDA
%       b(5).name = 'qda';   b(5).options.p = [];                             % QDA
%       b(6).name = 'nnglm'; b(6).options.method = 3; b(6).options.iter = 10; % Nueral network
%       b(7).name = 'svm';   b(7).options.kernel = 4;                         % rbf-SVM
%       b(8).name = 'maha';  b(8).options = [];                               % Euclidean distance
%       b(9).name = 'dmin';  b(9).options = [];                               % Mahalanobis distance
%       op = b;
%       op = Bcl_structure(X,d,op);                                           % Training only
%
%    Example: Testing only (after training only example):
%       ds = Bcl_structure(Xt,op);                                            % Testing only
%       p  = Bev_performance(ds,dt)
%
%    See also Bcl_exe.
%
% D.Mery, PUC-DCC, Jul 2009
% http://dmery.ing.puc.cl


function [ds,options] = Bcl_structure(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});

b = options;
n = length(b); % number of classifiers

if train
    for i=1:n
        Bn = b(i).name;
        if Bn(1) ~= 'B'
            Bn = ['Bcl_' Bn];
        end
        b(i).options = feval(Bn,X,d,b(i).options);
    end
    options = b;
    ds = options;
end
if test
    nt = size(Xt,1);
    ds3 = zeros(nt,n,2);
    d3 = 0;
    for i=1:n
        Bn = b(i).name;
        if Bn(1) ~= 'B'
            Bn = ['Bcl_' Bn];
        end
        dsi = feval(Bn,Xt,b(i).options);
        ds3(:,i,1) = dsi(:,1);
        if size(dsi,2)==2
            ds3(:,i,2) = dsi(:,2);
            d3 = 1;
        end
    end
    if d3
        ds = ds3;
    else
        ds = ds3(:,:,1);
    end
end

