% ds      = Bcl_lda(X,d,Xt,[])  Training & Testing together
% options = Bcl_lda(X,d,[])     Training only
% ds      = Bcl_lda(Xt,options) Testing only
%
% Toolbox: Balu
%    LDA (linear discriminant analysis) classifier.
%    We assume that the classes have a common covariance matrix
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.p is the prior probability, if p is empty,
%       it will be estimated proportional to the number of samples of each
%       class.
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.dmin contains min(d).
%       options.Cw1 is inv(within-class covariance).
%       options.mc contains the centroids of each class.
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'lda     ').
%
%    Reference:
%       Hastie, T.; Tibshirani, R.; Friedman, J. (2001): The Elements of
%       Statistical Learning, Springer (pages 84-90)
%
%    Example: Training & Test together:
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.p = [];
%       ds = Bcl_lda(X,d,Xt,op);   % LDA classifier
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.p = [0.75 0.25];        % prior probability for each class
%       op = Bcl_lda(X,d,op);      % LDA - training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_lda(Xt,op);       % LDA - testing
%       p = Bev_performance(ds,dt) % performance on test data
%
%    See also Bcl_qda.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl



function [ds,options] = Bcl_lda(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});

options.string = 'lda     ';
if train
    dmin = min(d);
    dmax = max(d);
    d    = d-dmin+1;
    N = length(d);   % number of samples
    K = dmax-dmin+1; % number of classes

    pest = isempty(options.p);

    if pest
        p = zeros(K,1);
    else
        p = options.p;
    end

    m   = size(X,2);
    L   = zeros(K,1);
    Cw = zeros(m,m);

    mc = zeros(m,K);
    for k=1:K
        ii   = find(d==k); % index of rows of class k
        if isempty(ii)
            error('Bcl_lda: There is no class %d in the data.',k+dmin-1)
        end
        L(k)      = length(ii); % number of samples in class k
        Xk        = X(ii,:);    % samples of class k
        mc(:,k)  = mean(Xk)';   % mean of class k
        Ck        = cov(Xk);    % covariance of class k
        Cw = Cw + Ck*(L(k)-1);  % within-class covariance
        if pest
            p(k) = L(k)/N;
        end
    end

    Cw = Cw/(N-K);
    options.Cw1  = inv(Cw);
    options.dmin = dmin;
    options.mc   = mc;
    options.p    = p;
    ds = options;
end

if test
    K  = size(options.mc,2);
    Nt = size(Xt,1);
    D  = zeros(Nt,K);
    for k=1:K
        C1 = options.Cw1*options.mc(:,k);
        C2 = (-0.5*options.mc(:,k)'*C1 + log(options.p(k)))*ones(Nt,1);
        D(:,k) = Xt*C1+C2;
    end
    [i,j]=max(D,[],2);
    sc = ones(size(i))./(abs(i)+1e-5);
    ds = j + options.dmin - 1;
    ds = Bcl_outscore(ds,sc,options);
end

