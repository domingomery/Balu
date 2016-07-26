% ds      = Bcl_dmin(X,d,Xt,[])  Training & Testing together
% options = Bcl_dmin(X,d,[])     Training only
% ds      = Bcl_dmin(Xt,options) Testing only
%
% Toolbox: Balu
%    Classifier using Euclidean minimal distance
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
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'dmin    ').
%
%    Example: Training & Test together:
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       ds = Bcl_dmin(X,d,Xt,[]);  % Euclidean distance classifier
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op = Bcl_dmin(X,d,[]);     % Euclidean distance classifier - training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_dmin(Xt,op);      % Euclidean distance classifier - testing
%       p = Bev_performance(ds,dt) % performance on test data
%
%    See also Bcl_maha.
%
% D.Mery, PUC-DCC, May 2010
% http://dmery.ing.puc.cl


function ds = Bcl_dmin(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});

options.string = 'dmin    ';
if train
    m    = size(X,2);
    dmin = min(d);
    dmax = max(d);
    d    = d-dmin+1;
    n    = dmax-dmin+1;
    mc   = zeros(n,m);
    for i=1:n
        mc(i,:) = mean(X(d==i,:),1);
    end
    options.mc   = mc;
    options.dmin = dmin;
    ds = options;
end
if test
    mc = options.mc;
    n  = size(mc,1);
    Nt = size(Xt,1);
    ds = zeros(Nt,1);
    sc = zeros(Nt,1);
    for q=1:Nt
        D     = ones(n,1)*Xt(q,:)-mc;
        D2    = D.*D;
        e     = sum(D2,2);
        [i,j] = min(e);
        ds(q) = j;
        sc(q) = i;
    end
    ds = ds+options.dmin-1;
    ds = Bcl_outscore(ds,sc,options);
end


