% s = Bex_sfssvm(f,d,m)
%
% Toolbox: Balu
%    Example: Feature selection using Balu algorithm based on SFS and SVM.
%
%    Bfs_balu has three steps:
%       (1) normalizes (using Bft_norm), 
%       (2) cleans (using Bfs_clean), and 
%       (3) selects features (using Bfs_sfs).
%
% Example:
%   load datareal
%   s = Bex_sfssvm(f,d,6)
%   fn(s,:)
%
% See also Bfs_balu
% 
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl


function s = Bex_sfssvm(f,d,m)

op.m = m;                      % m features will be selected
op.s = 0.75;                   % only 75% of sample will be used
op.show = 1;                   % display results
op.b.name = 'svm';             % SFS with SVM
op.b.options.kernel = 4;       % SVM-RBF
s = Bfs_balu(f,d,op);          % index of selected features
X = f(:,s);                    % selected features
op = op.b.options;
[X1,d1,X2,d2] = Bds_stratify(X,d,0.75);
d2s = Bcl_svm(X1,d1,X2,op);
p = Bev_performance(d2,d2s);     % performance with lsef
fprintf('Performance with SFS and SVM = %5.4f\n',p)
