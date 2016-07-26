% ds      = Bcl_maha(X,d,Xt,[])  Training & Testing together
% options = Bcl_maha(X,d,[])     Training only
% ds      = Bcl_maha(Xt,options) Testing only
%
% Toolbox: Balu
%    Classifier using Mahalanobis minimal distance
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.mc contains the centroids of each class.
%       options.dmin contains min(d).
%       options.Ck is covariance matrix of each class.
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'maha    ').
%
%    Example: Training & Test together:
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       ds = Bcl_maha(X,d,Xt,[]);  % Euclidean distance classifier
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op = Bcl_maha(X,d,[]);     % Euclidean distance classifier - training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_maha(Xt,op);      % Euclidean distance classifier - testing
%       p = Bev_performance(ds,dt) % performance on test data
%
%    See also Xdmin.
%
% D.Mery, PUC-DCC, May 2010
% http://dmery.ing.puc.cl

function [ds,options] = Bcl_maha(varargin)
[train,test,X,d,Xt,options] = Xconstruct(varargin{:});

options.string = 'maha    ';
if train
    m    = size(X,2);
    dmin = min(d);
    dmax = max(d);
    n    = dmax-dmin+1;
    d    = d-dmin+1;
    mc   = zeros(n,m);
    M    = size(X,2);
    Ck   = zeros(M,M,n);
    for i=1:n
        ii = find(d==i);
        mc(i,:) = mean(X(ii,:),1);
        CCk = cov(X(ii,:));      % covariance of class i
        Ck(:,:,i) = CCk;
    end
    options.mc   = mc;
    options.dmin = dmin;
    options.Ck    = Ck;
    ds = options;
end
if test
    mc    = options.mc;
    n     = size(mc,1);
    Nt    = size(Xt,1);
    ds    = zeros(Nt,1);
    sc    = ds;
    for k=1:n
        Ck(:,:,k) = inv(options.Ck(:,:,k));
    end
    for q=1:Nt
        dk = zeros(n,1);
        for k=1:n
            dx = Xt(q,:)-options.mc(k,:);
            dk(k) = dx*Ck(:,:,k)*dx';
        end
        [i,j] = min(dk);
        ds(q) = j;
        sc(q) = i;
    end
    ds = ds+options.dmin-1;
    ds = Xoutscore(ds,sc,options);
end

