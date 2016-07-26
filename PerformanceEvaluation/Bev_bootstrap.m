% [p,ci] = Bev_bootstrap(X,d,options)
%
% Toolbox: Balu
%    Bootstrap evaluation in B bootstrap samples of X and classification d
%    according to given classifier.
%
%    X is a matrix with features (columns)
%    d is the ideal classification for X
%
%    options.b is a Balu classifier or several classifiers (see example)
%    options.B is the number bootstrap samples
%    options.c is the probability of the confidence intervale.
%    options.show = 1 displays results.
%
%    Example for one classifier:
%       load datagauss                                      % simulated data (2 classes, 2 features)
%       Bplotfeatures(X,d)                                  % plot feature space
%       b.name = 'knn'; b.options.k = 5;                    % knn with 5 neighbors
%       op.b = b; op.B = 50; op.show = 0; op.c = 0.90;      % evaluation in 30 bootstrap samples
%       [p,ci] = Bev_bootstrap(X,d,op)                      % Bootstrap
%
%    Example for more classifiers:
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
%       op.b = b; op.B = 30; op.show = 1; op.c = 0.95;                        % evaluation in 30 bootstrap samples
%       [p,ci] = Bev_bootstrap(X,d,op);                                       % Bootstrap
%
%    The mean performance of classifier k is given in p(k), and the
%    confidence intervals for c*100% are in ci(k,:).
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [p,ci] = Bev_bootstrap(X,d,options)

B    = options.B;
b    = options.b;
c    = options.c;
show = options.show;
Xr   = X;
dr   = d;
dmin = min(d);
dmax = max(d);
nn   = dmin:dmax;
n    = length(b);
th   = zeros(B,n);
p    = zeros(n,1);
ci   = zeros(n,2);

for i=1:B    
    [XXb,ddb,XXnb,ddnb] = Bds_bootstrap(Xr,dr);
    [dds,ops]           = Bcl_structure(XXb,ddb,XXnb,b);
    th(i,:)             = Bev_performance(ddnb,dds,nn);
end

for k=1:n
    s = ops(k).options.string;
    thk = th(:,k);
    [i,j] = sort(thk);
    c1 = fix(B*(1-c)/2+0.5);
    c2 = B-c1;
    p(k) = mean(thk);
    p1   = thk(j(c1));
    p2   = thk(j(c2));
    ci(k,:) = [p1 p2];
    if show
        fprintf('%3d) %s  %5.2f%% in (%5.2f, %5.2f%%) with CI=%2.0f%% \n',k,s,p(k)*100,p1*100,p2*100,c*100);
    end
end