% s = Bex_exsknn(f,d,m)
%
% Toolbox: Balu
%    Example: Feature selection using exhaustive search and KNN.
%
%    Bfs_balu has three steps:
%       (1) normalizes (using Bft_norm),
%       (2) cleans (using Bfs_clean), and
%       (3) selects features (using Bfs_sfs).
%
% Example:
%   load datareal
%   s = Bex_exsknn(f,d)
%   fn(s,:)
%
% See also Bfs_balu
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl


function s = Bex_exsknn(f,d)

op.m = 10;                       % 10 features will be pre-selected
op.s = 0.75;                     % only 75% of sample will be used
op.show = 1;                     % display results
op.b.name = 'fisher';            % SFS with fisher
s1 = Bfs_balu(f,d,op);           % index of selected features
X1 = f(:,s1);                    % selected features
op.m = 4;                        % 3 features will be selected
op.show = 1;                     % display results
op.b.name = 'knn';               % SFS with KNN
op.b.options.k = 5;              % 5 neighbors
s2 = Bfs_exsearch(X1,d,op);      % index of selected features
s = s1(s2);                      % list of feature names
X = f(:,s);
op = op.b.options;
[X1,d1,X2,d2] = Bds_stratify(X,d,0.75);
d2s = Bcl_knn(X1,d1,X2,op);
p = Bev_performance(d2,d2s);     % performance with lsef
fprintf('Performance with exhaustive search and KNN = %5.4f\n',p)
