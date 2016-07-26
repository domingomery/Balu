% [T,p] = Bev_holdout(X,d,options)
%
% Toolbox: Balu
%    Holdout evaluation of a classifier.
%
%    X is a matrix with features (columns)
%    d is the ideal classification for X
%
%    options.b is a Balu classifier or several classifiers (see example)
%    options.s is the portion of data used for training, e.g. s=0.75.
%    options.strat = 1 means the portion is stratified.
%
%    Confusion matrix is given file by file for each group in matrix T.
%    The performance in each group is given vector per.
%
%    Example for one classifier:
%       load datagauss                      % simulated data (2 classes, 2 features)
%       Bplotfeatures(X,d)                  % plot feature space
%       b.name = 'knn'; b.options.k = 5;   % knn with 5 neighbors
%       op.b = b; op.strat = 1;op.s = 0.75; % stratify with 75% train
%       [T,p] = Bev_holdout(X,d,op)         % holdout
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
%       op.b = b; op.strat = 1; op.s = 0.75;                                  % stratify with 75% train
%       [T,p] = Bev_holdout(X,d,op);                                          % holdout
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [T,p] = Bev_reclassification(X,d,options)

if not(exist('show','var'))
    show=1;
end

b = options.b;

[ds,op] = Bcl_structure(X,d,X,b);

nn   = min(d):max(d);
m = length(nn);
n    = length(b);
T = zeros(m,m,n);
p = zeros(n,1);

for i=1:n
    [TT,pp] = Bev_confusion(d,ds(:,i),nn);
    T(:,:,i) = TT;
    p(i) = pp;
    if show
        s = op(i).options.string;
        fprintf('%3d) %s  %5.2f%% \n',i,s,p(i)*100);
    end
end

