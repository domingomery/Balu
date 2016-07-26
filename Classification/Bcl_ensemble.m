% ds      = Bcl_ensemble(X,d,Xt,options)  Training & Testing together
% options = Bcl_ensemble(X,d,options)     Training only
% ds      = Bcl_ensemble(Xt,options)      Testing only
%
% Toolbox: Balu
%    Design and test an ensemble of n classifiers.
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       tensemble is the type of the ensemble:
%          'vote'  for majority vote
%          'min'   for minimum
%          'max'   for maximum
%          'best'  takes the best individual classifier reclassifying X
%          'all'   for all classifications (one classification per column in ds)
%          'stack' uses a stacked generalization acoording to Polikar-Fig 11.
%                  in this case param is the classifier defined as b
%                  structure.
%
%       b is the following structure for classifier i (i=1,...,n)
%          b(i).name = classifier's name
%          b(i).p1   = parameter 1 of the classifier if any
%          b(i).p2   = parameter 2 of the classifier if any
%          b(i).p3   = parameter 3 of the classifier if any
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%
%    show = 1 displays results (default show = 0).
%
%    Polikar, R. (2006): Ensemble Based Systems in Decision Making. IEEE
%    Circuits ans Systems Magazine, Third Quarter, 21-45.
%
%    Giacinto, G. & Roli, F. (2001): Dynamic classifier selection based on
%    multiple classifier behaviour Pattern Recognition, 34: 1879 - 1881
%
%    Example: Training & Test together:
%       load datagauss
%       b(1).name = 'knn';   b(1).options.k = 5;                              % KNN with 5 neighbors
%       b(2).name = 'knn';   b(2).options.k = 7;                              % KNN with 7 neighbors
%       b(3).name = 'knn';   b(3).options.k = 9;                              % KNN with 9 neighbors
%       b(4).name = 'lda';   b(4).options.p = [];                             % LDA
%       b(5).name = 'qda';   b(5).options.p = [];                             % QDA
%       b(6).name = 'nnglm'; b(6).options.method = 3; b(6).options.iter = 10; % Nueral network
%       b(7).name = 'svm';   b(7).options.kernel = 4;                         % rbf-SVM
%       b(8).name = 'maha';  b(8).options = [];                               % Euclidean distance
%       b(9).name = 'dmin';  b(9).options = [];                               % Mahalanobis distance
%       op.b   = b;
%       op.tensemble = 'vote';
%       op.show = 1;
%       op.bs   = [];
%       ds1 = Bcl_ensemble(X,d,Xt,op);               % majority vote
%       p1  = Bev_performance(ds1,dt)                % performance on test data
%
%       op.tensemble = 'max';
%       op.bs   = [7 0.8 0.7];
%       ds2 = Bcl_ensemble(X,d,Xt,op);               % DCS
%       p2  = Bev_performance(ds2,dt)                % performance on test data
%
%
%       op.tensemble = 'stack';                      % STACK
%       bs.name = 'svm'; bs.options.kernel = 4;      % rbf-SVM
%       op.bs = bs;
%       ds3 = Bcl_ensemble(X,d,Xt,op);               % STACK
%       p3  = Bev_performance(ds3,dt)                % performance on test data
%
%    Example: Training only
%       load datagauss
%       b(1).name = 'knn';   b(1).options.k = 5;                              % KNN with 5 neighbors
%       b(2).name = 'knn';   b(2).options.k = 7;                              % KNN with 7 neighbors
%       b(3).name = 'knn';   b(3).options.k = 9;                              % KNN with 9 neighbors
%       b(4).name = 'lda';   b(4).options.p = [];                             % LDA
%       b(5).name = 'qda';   b(5).options.p = [];                             % QDA
%       b(6).name = 'nnglm'; b(6).options.method = 3; b(6).options.iter = 10; % Nueral network
%       b(7).name = 'svm';   b(7).options.kernel = 4;                         % rbf-SVM
%       b(8).name = 'maha';  b(8).options = [];                               % Euclidean distance
%       b(9).name = 'dmin';  b(9).options = [];                               % Mahalanobis distance
%       op.b   = b;
%       op.tensemble = 'vote';
%       op.show = 1;
%       op.bs   = [];
%       op  = Bcl_ensemble(X,d,op);               % training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_ensemble(Xt,op);                 % testing
%       p = Bev_performance(ds,dt) 
%
%   See also Bcl_structure, Bcl_dcs.
%
% D.Mery, PUC-DCC, Jul 2009
% http://dmery.ing.puc.cl
%


function [ds,options] = Bcl_ensemble(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'ensemble';
Bsclass   = options.b;
tensemble = options.tensemble;
bs        = options.bs;
show      = options.show;
if train
    opt = Bcl_structure(X,d,Bsclass);
    switch lower(tensemble)
        case 'best'
            DS = Bcl_structure(X,d,X,Bsclass);
            p = Bev_performance(d,DS);
            [i,j] = max(p);
            options.j = j;
        case 'stack'
            DS = Bcl_structure(X,d,X,Bsclass);
            op = bs;
            ops = Bcl_structure(DS,d,op);
            options.ops = ops;
    end
    options.opt = opt;
    ds = options;
end
if test
    opt = options.opt;
    DST = Bcl_structure(Xt,opt);
    switch lower(tensemble)
        case 'max'
            ds = max(DST,[],2);
        case 'min'
            ds = min(DST,[],2);
        case 'vote'
            ds = mode(double(DST),2);
        case 'all'
            ds = DST;
        case 'best'
            j = options.j;
            ds = DST(:,j);
            if show
                fprintf('Best classifier is number %d (%s) with performance = %5.2f%%.n',j,Bsclass(j).name,p(j)*100');
            end
        case 'stack'
            ops = options.ops;
            ds = Bcl_structure(DST,ops);
%       'dcs'   uses dynamic classifier selection based on multiple
%               classifier behaviour after Giacinto (2001).
%        case 'dcs'
%            op.b = Bsclass;
%            op.k    = bs(1);
%            op.ths  = bs(2);
%            op.thc  = bs(3);
%            ds = Bdcs(X,d,Xt,op);
        otherwise
            error('Bcl_structure does not recognize %s type of ensemble.',tensemble)
    end

end






