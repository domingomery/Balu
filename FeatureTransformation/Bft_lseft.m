% [Y,selec,th] = Bft_lseft(X,d,options)
%
% Toolbox: Balu
%    Feature transformacion using LSE-forward algorithm
%
%    input: X feature matrix
%           options.m number of features to be selected
%           optoins.show = 1 displays results
%           options.pca = 1 for PCA and = 0 for PLS
%
%    output: selec selected features
%            Y is equal to A*th, where A = [X(:,selec) ones(size(X,1),1)]
%            Y is very similar to PCA or PLS.
%
%    The main idea of the algorithms is to evaluate and select feature
%    subsets based on their capacities to reproduce sample projections on
%    principal axes:
%
%    Reference:
%    Mao, K. Identifying critical variables of principal components for
%    unsupervised feature selection Systems, Man, and Cybernetics, Part B:
%    Cybernetics, IEEE Transactions on, 2005, 35, 339-344
%
%    Example using PLSR
%    load datareal
%    op.m = 20;                      % 20 features will be pre-selected
%    op.show = 1;                    % display results
%    op.s = 1;
%    op.b.name = 'fisher';           % SFS with Fisher
%    s0 = Bfs_balu(f,d,op);          % index of selected features
%    X0 = f(:,s0);                   % preselected features
%    op.m = 10;                      % 10 features will be selected
%    op.show = 1;                    % display results
%    X1 = Bft_plsr(X0,d,op.m);
%    op.pca = 0;
%    X2 = Bft_lseft(X0,d,op);        % index of selected features
%    X3 = X0(:,1:op.m);
%    op.p = [];
%    ds1 = Bcl_lda(X1,d,X1,op);
%    p1 = Bev_performance(d,ds1);    % performance with PLSR
%    ds2 = Bcl_lda(X2,d,X2,op);
%    p2 = Bev_performance(d,ds2);    % performance with LSEFT
%    ds3 = Bcl_lda(X3,d,X3,op);
%    p3 = Bev_performance(d,ds3);    % performance with SFS
%    fprintf('Performance with %d features: SFS = %5.4f, with PLSR = %5.4f, with LSEFT = %5.4f\n',op.m,p3,p1,p2)
%
% See also Bfs_lsef
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [Y,selec,th] = Bft_lseft(X,d,options)

[selec,Y,th] = Bfs_lsef(X,d,options);


