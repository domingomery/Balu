% selec = Bfs_mRMR(X,d,options)
%
% Toolbox: Balu
%    Feature selection using Criteria of Max-Dependency, Max-Relevance, and
%    Min-Redundancy after Peng et al. (2005)
%
%    X extracted features (NxM): N samples, M features
%    d ideal classification (Nx!) for the N samples
%    options.m number of selected features
%    options.p a priori probability of each class (if not given, is
%    estimated using the ratio samples per class to N.
%
% References:
%   Peng, H.; Long, F.; Ding, C. (2005): Feature Selection Based on Mutual 
%   Information: Criteria of Max-Dependency, Max-Relevance, and Min-
%   Redundancy. IEEE Trans. on Pattern Analysis and Machine Intelligence, 
%   27(8):1226-1238.
%
%   Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010): Kernel density 
%   estimation via diffusion Annals of Statistics, 38(5):2916-2957. 
%
% NOTE:
%    The pdf's are estimated using Kernel Density Estimations programs 
%    kde.m and kde2d.m after Botev et al. (2010) implemented by Botev. 
%    These files are in Balu directory 'Feature Analysis' as Bfa_kde and
%    Bfs_kde2d. They can also be downloaded from www.mathwork.com
%    (c) Zdravko Botev. All rights reserved.
%
% Example (comparison between SFS-Fisher and mRMR): 
%
%    load datareal
%    s = Bfs_clean(f,1);
%    X = f(:,s);
%    k = 0;
%    k=k+1;b(k).name = 'knn';     b(k).options.k = 9;               % KNN with 9 neighbors
%    k=k+1;b(k).name = 'lda';     b(k).options.p = [];              % LDA
%    k=k+1;b(k).name = 'qda';     b(k).options.p = [];              % QDA
%    k=k+1;b(k).name = 'dmin';    b(k).options = [];                % Euclidean distance
%    k=k+1;b(k).name = 'maha';    b(k).options = [];                % Mahalanobis distance
%    k=k+1;b(k).name = 'libsvm';  b(k).options.kernel = '-t 2';                    % rbf-SVM
%    k=k+1;b(k).name = 'nnglm';   b(k).options.method = 3; b(k).options.iter = 10; % Nueral network
%
%    op.strat=1; op.b = b; op.v = 10; op.show = 1; op.c = 0.95;     % 10 groups cross-validation
%
%    msel1 = 12;
%    msel2 = 4;
%    
%    op1.m = msel1;op1.show=1;op1.b.name='fisher';
%    s1    = Bfs_sfs(X,d,op1);
%    X1    = X(:,s1);
%    disp('Performances using SFS-Fisher:')
%    p1    = Bev_crossval(X1(:,1:msel2),d,op);
%    
%    op2.m = msel2;op2.show=1;                    
%    s2    = Bfs_mRMR(X1,d,op2);
%    disp(' ')
%    disp('Performances using SFS-mRMR:')
%    p2 = Bev_crossval(X1(:,s2),d,op); 
%    disp(' ')
%    disp('Comparison of performances and mean of performances:')
%    [p1 p2], mean([p1 p2])
%
% D.Mery, PUC-DCC, Jul. 2011
% http://dmery.ing.puc.cl

function s = Bfs_mRMR(X,d,options)

if ~isfield(options,'p')
    dn = max(d)-min(d)+1; % number of classes
    p = ones(dn,1)/dn;
else
    p = options.p;
end

if ~isfield(options,'show')
    show = 0;
else
    show = options.show;
end

% T = 100;

M  = size(X,2);
n  = options.p;
s  = zeros(n,1);
t  = ones(M,1);
h  = zeros(n,1);
ff = Bio_statusbar('Bfs_mRMR');

Ic = zeros(M,1);
m  = 1;
for j=1:M
    %   Ic(j) = Bfa_mutualinfo(X(:,j),d,p); % see in paper I(xj;c)
    Ic(j) = Bfa_miparzen2([X(:,j) d],1,p); % see in paper I(xj;c)
end
ff = Bio_statusbar(1/n,ff);
[Icmax,jsel] = max(Ic);
if show
    clf
    h(m) = Icmax;
    bar(h);
end

s(m)    = jsel;
t(jsel) = 0;


MI = zeros(M,M);

if n>1
    for m=2:n
        
        Jmax = -Inf;
        for j=1:M
            if t(j)==1 % if yes, feature j is not selected 
                xj = X(:,j);
                sumI = 0;
                for i=1:m-1
                    xi = X(:,s(i));
                    % sumI = sumI + Bfa_mutualinfo2([xj xi]); % see in paper I(xj;xi)
                    if MI(s(i),j)==0;
                        mij = Bfa_miparzen2([xj xi]);
                        MI(s(i),j) = mij;
                        MI(j,s(i)) = mij;
                    else
                        mij = MI(s(i),j);
                    end
                    sumI = sumI + mij; % see in paper I(xj;xi)
%                    sumI = sumI + Bfa_miparzen2([xj xi]); % see in paper I(xj;xi)
                end
                % Jj = Ic(j)-sumI/(m-1);
                Jj = Ic(j)/(sumI/(m-1));
                % Jj = Ic(j)-sumI;
                if Jj>Jmax
                    Jmax = Jj;
                    jsel = j;
                end
            end
        end
        s(m)    = jsel;
        t(jsel) = 0;
        if show
            h(m) = Jmax;
            bar(h)
        end
    ff = Bio_statusbar(m/n,ff);
    end
end
delete(ff)


