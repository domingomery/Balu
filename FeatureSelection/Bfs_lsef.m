% [selec,Y,th] = Bfs_lsef(X,options)
%
% Toolbox: Balu
%    Feature Selection using LSE-forward algorithm
%
%    input: X feature matrix
%           options.m number of features to be selected
%           optoins.show = 1 displays results
%
%    output: selec selected features
%            Y is equal to A*th, where A = [X(:,selec) ones(size(X,1),1)]
%            Y is very similar to PCA.
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
%    Example 1: using PCA
%    load datareal
%    s1 = [279 235 268 230 175 165 207 160 269 157]; %indices using Example
%                                                    %of Bfs_sfs.
%    X1 = f(:,s1);                   % preselected features
%    Xn1 = fn(s1,:);
%    op.m = 6;                       % 5 features will be selected
%    op.show = 1;                    % display results
%    op.pca = 1;
%    s = Bfs_lsef(X1,d,op);          % index of selected features
%    T1 = X1(:,s);                   % selected features (transformation
%                                    %                     is avoided)
%    op.p = [];
%    ds1 = Bcl_lda(T1,d,T1,op);
%    p1 = Bev_performance(d,ds1)     % performance with lsef
%    T2 = Bft_pca(X1,op.m);          % transformed features using PCA
%    op.p = [];
%    ds2 = Bcl_lda(T2,d,T2,op);
%    p2 = Bev_performance(d,ds2)     % performance with PCA
%    fprintf('Performance with lsef = %5.4f, with PCA = %5.4f\n',p1,p2)
%
%    Example 2: Linear transformation of selected features is similar
%               to PCA.
%    load datareal
%    s1 = [279 235 268 230 175 165 207 160 269 157]; %indices using Example
%                                                    %of Bfs_sfs.
%    X1 = f(:,s1);                   % preselected features
%    Xn1 = fn(s1,:);
%    op.m = 6;                       % 6 features will be selected
%    op.show = 1;                    % display results
%    op.pca = 1;
%    [s,Y,th] = Bfs_lsef(X1,d,op);   % Y is computed as A*th, where
%                                    % A = [X1(:,s) ones(size(X1,1),1)]
%    Yt = Bft_pca(X1,op.m);          % PCA analysis
%    sqdif(Y,Yt)                       % Y and Yt are very similar.
%
% See also Bft_lseft, Bft_pca, Bft_plsr.
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [selec,Y,th] = Bfs_lsef(X,d,options)

% d is not used by this algorithm.
% the sintaxis "Bfs_fosmod(X,d,options)" is allowed in order to be
% similar to other Bfs_ functions.

if nargin==2;
    options = d;
end

m    = options.m;
show = options.show;
if isfield(options,'pca')
    pcat = options.pca;
else
    pcat = 0;
end


[N,n] = size(X); % number of instances (samples), number of features

% 1) Initialize S to an empty set
% S = [];

% 2) Initialize R to the full feature set
% R = X;

% 3) Perform PCA on complete data:
% 4) Calaculate sample projections on the first principal axes with 0.95
%    energy
if pcat
    Yt = Bft_pca(X,m);
else
    Yt = Bft_plsr(X,d,m);
end
Ny   = length(Yt(:));

% 5) m iterations for m features
p  = ones(n,1); % '1' means not selected feature
for l=1:m
    em = Inf*ones(n,1);
    S  = X(:,not(p));
    Yss = zeros(N,Ny/N,n);
    ths = zeros(l+1,Ny/N,n);
    for j=1:n
        if p(j)
            A     = [S X(:,j) ones(N,1)];
            th    = (A'*A)\A'*Yt;
            ths(:,:,j) = th;
            Ys    = A*th;
            em(j) = norm(Yt-Ys);
            Yss(:,:,j) = Ys;
        end
    end
    [emin,js] = min(em); % js is the feature that minimizes Yt-Ys, i.e.,
    % Ys is the best fit of Yt computed as a linear
    % transformation of [S X(:,js) ones(N,1)]
    Y = Yss(:,:,js);
    th = ths(:,:,js);
    if show
        fprintf('%2d) selected feature=%4d error=%7.2f%%\n',l,js,emin/Ny*100)
    end
    p(js) = 0;
end
selec = find(p==0);
