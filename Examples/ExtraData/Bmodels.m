function [bfs,bcl] = Bmodels()

% Feature Selection Models
bfs(1).name = 'sfs';   bfs(1).options.b.name    = 'fisher';
bfs(2).name = 'sfs';   bfs(2).options.b.name    = 'knn'; bfs(2).options.b.options.k = 5;
bfs(3).name = 'sfs';   bfs(3).options.b.name    = 'sp100';
bfs(4).name = 'sfs';   bfs(4).options.b.name    = 'knn'; bfs(4).options.b.options.k = 7;
bfs(5).name = 'rank';  bfs(5).options.criterion = 'roc';

% Classifier Models
bcl(1).name = 'knn';   bcl(1).options.k = 5;      % KNN with 5 neighbors
bcl(2).name = 'knn';   bcl(2).options.k = 7;      % KNN with 7 neighbors
bcl(3).name = 'knn';   bcl(3).options.k = 9;      % KNN with 9 neighbors
bcl(4).name = 'lda';   bcl(4).options.p = [];     % LDA
bcl(5).name = 'qda';   bcl(5).options.p = [];     % QDA
bcl(6).name = 'svm';   bcl(6).options.kernel = 4; % rbf-SVM
bcl(7).name = 'maha';  bcl(7).options = [];       % Euclidean distance
bcl(8).name = 'dmin';  bcl(8).options = [];       % Mahalanobis distance
