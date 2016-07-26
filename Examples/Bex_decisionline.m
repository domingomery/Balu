% Bex_decisionline(X,d,bcl)
%
% Toolbox: Balu
%    Example: Decision lines of features X with labels d for classifiers 
%    defined in bcl
%
% Example:
%    load datagauss                                    % simulated data (2 classes, 2 features)
%    Xn = ['x_1';'x_2'];
%    bcl(1).name = 'knn';   bcl(1).options.k = 5;          % KNN with 5 neighbors
%    bcl(2).name = 'lda';   bcl(2).options.p = [];         % LDA
%    bcl(3).name = 'svm';   bcl(3).options.kernel = 3;     % rbf-SVM
%    Bex_decisionline(X,d,bcl)
%
% See also Bio_decisionline
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl

function Bex_decisionline(X,d,bcl)

op = Bcl_structure(X,d,bcl);
Bio_decisionline(X,d,['x1';'x2'],op);
