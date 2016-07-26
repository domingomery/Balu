% [T,p] = Bev_jackknife(X,d,options)
%
% Toolbox: Balu
%    Holdout evaluation of a classifier.
%
%    v-fold Cross Validation in v groups of samples X, where v is the
%    number of sanples (i.e., v = size(X,1)). The training will be in X 
%    without sample i and testing in sample i. ci is the confidence interval.
%
%    X is a matrix with features (columns)
%    d is the ideal classification for X
%
%    options.b is a Balu classifier or several classifiers (see example)
%    options.c is the probability of the confidence intervale.
%    options.show displays results.
%
%    Example for one classifier:
%       load datagauss                      % simulated data (2 classes, 2 features)
%       Bplotfeatures(X,d)                  % plot feature space
%       b.name = 'knn'; b.options.k = 5;    % knn with 5 neighbors
%       op.b = b; op.c = 0.90; op.show = 0; % confidence intervale for 90
%       [p,ci] = Bev_jackknife(X,d,op)      % jackknife
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
%       b(8).name = 'maha';  b(8).options = [];                               % Mahalanobis distance
%       b(9).name = 'dmin';  b(9).options = [];                               % Euclidean distance
%       op.b = b; op.show = 1; op.c = 0.95;                                   % 95% confidence interval
%       [p,ci] = Bev_jackknife(X,d,op);                                       % jackknife
%
%    The mean performance of classifier k is given in p(k), and the
%    confidence intervals for c*100% are in ci(k,:).
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [T,p] = Bev_jackknife(X,d,options)
options.v = size(X,1);
options.strat = 0;
[T,p] = Bev_crossval(X,d,options);
