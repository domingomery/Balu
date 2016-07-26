% ds      = Bcl_bagging(X,d,Xt,options)  Training & Testing together
% options = Bcl_bagging(X,d,options)     Training only
% ds      = Bcl_bagging(Xt,options)      Testing only
%
% Toolbox: Balu
%    Bagging classifier.
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.b is a Balu classifier
%       options.Nboost is the number of boosting classifiers
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.dmin contains min(d).
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'bagging ').
%       options.opt is an array with options.NBoost obtions for each
%       classifier specified by options.b.
%
%    See reference:
%    Witten, I.H.; Frank, E. (2005): Data Mining, Elsevier, 2nd Edition
%
%    Example:
%       load datagauss               % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)        % plot feature space
%       b.name = 'knn';
%       b.options.k = 5;             % definition of KNN with k=5
%       op.b = b;
%       op.NBoost = 15;              % Bagging with 15 boostings
%       ds = Bcl_bagging(X,d,Xt,op); % Bagging with knn (5 neighbors)
%       p = Bev_performance(ds,dt)   % performance on test data
%
%    Example: Training & Test together:
%       load datagauss               % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)        % plot feature space
%       b.name = 'knn';
%       b.options.k = 5;             % definition of KNN with k=5
%       op.b = b;
%       op.NBoost = 15;              % Bagging with 15 boostings
%       op = Bcl_bagging(X,d,op);    % Bagging with knn (5 neighbors)
%
%    Example: Training only
%       ds = Bcl_bagging(Xt,op);     % Bagging with knn (5 neighbors)
%       p = Bev_performance(ds,dt)   % performance on test data
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl



function [ds,options] = Bcl_bagging(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'bagging ';
if train
    b = options.b;
    B = options.NBoost;
    % op.b = b;
    N = size(X,1);
    dmin = min(d);
    d = d-dmin+1;
    for i=1:B
        i1  = ceil(N*rand(N,1));
        XX = X(i1,:);
        dd = d(i1);
        opt(i) = Bcl_structure(XX,dd,b); %#ok<AGROW>
    end
    options.dmin = dmin;
    options.opt = opt;
    ds = options;
end
if test
    nt   = size(Xt,1);
    ds   = zeros(nt,1);
    B    = options.NBoost;
    dmin = options.dmin;
    opt  = options.opt;
    dsb  = zeros(nt,B);
    for i=1:B
        dsb(:,i) = Bcl_structure(Xt,opt(i));
    end    
    for q=1:nt
        dv = dsb(q,:);
        ds(q)= mode(dv);
    end
    ds = ds+dmin-1;
end









