% s = Bex_lsef(f,d,m)
%
% Toolbox: Balu
%    Example: Feature selection using lsef algorithm
%
% Example:
%   load datareal
%   s = Bex_lsef(f,d,6)
%   fn(s,:)
%
% See also Bfs_lsef
% 
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl


function s = Bex_lsef(f,d,m)


op.m = 2*m;                     % 2*m features will be selected using SFS
op.s = 0.75;                    % only 75% of sample will be used
op.show = 1;                    % display results
op.b.name = 'fisher';           % SFS with Fisher
s1 = Bfs_balu(f,d,op);          % index of selected features
X1 = f(:,s1);                   % preselected features
op.m = m;                       % m/2 features will be selected
op.show = 1;                    % display results
[s2,Y,th] = Bfs_lsef(X1,op);    % Y is computed as A*th, where
                                % A = [X1(:,s2) ones(size(X1,1),1)]
Ypca = Bft_pca(X1,0.95);        % PCA analysis
mse(Y-Ypca)                     % Y and Ypca are very similar.

T1 = X1(:,s2);                  % selected features (transformation
op.p = [];
ds1 = Bcl_lda(T1,d,T1,op);
p1 = Bev_performance(d,ds1);    % performance with lsef
T2 = Bft_pca(X1,op.m);          % transformed features using PCA
op.p = [];
ds2 = Bcl_lda(T2,d,T2,op);
p2 = Bev_performance(d,ds2);    % performance with PCA
fprintf('Performance with lsef = %5.4f, with PCA = %5.4f\n',p1,p2)
s = s1(s2);
