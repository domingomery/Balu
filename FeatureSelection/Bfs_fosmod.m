% selec = Bfs_fosmod(X,options)
%
% Toolbox: Balu
%    Feature Selection using FOS-MOD algorithm
%
%    input: X feature matrix
%           options.m number of features to be selected
%           options.show = 1 displays results
%
%    output: selec selected features
%
%    FOS_MOD is a forward orthogonal search (FOS) algorithm by maximizing
%    the overall dependency (MOD). The algorithm assumes that a linear
%    relationship exists between sample features:
%
%    Reference:
%    Wei, H.-L. & Billings, S. Feature Subset Selection and Ranking for
%    Data Dimensionality Reduction Pattern Analysis and Machine
%    Intelligence, IEEE Transactions on, 2007, 29, 162-166
%
%    Example:
%    load datareal
%    s1 = [279 235 268 230 175 165 207 160 269 157]; %indices using Example
%                                                    %of Bfs_sfs.
%    X1 = f(:,s1);                   % preselected features
%    Xn1 = fn(s1,:);
%    op.m = 5;                       % 5 features will be selected
%    op.show = 1;                    % display results
%    s = Bfs_fosmod(X1,op);          % index of selected features
%    T1 = X1(:,s);                   % selected features (transformation
%                                    %                     is avoided)
%    op.p = [];
%    ds1 = Bcl_lda(T1,d,T1,op);
%    p1 = Bev_performance(d,ds1)     % performance with fosmod
%    T2 = Bft_pca(X1,op.m);          % transformed features using PCA
%    op.p = [];
%    ds2 = Bcl_lda(T2,d,T2,op);
%    p2 = Bev_performance(d,ds2)     % performance with PCA
%    fprintf('Performance with fosmod = %5.4f, with PCA = %5.4f\n',p1,p2)
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function selec = Bfs_fosmod(X,d,options)

% d is not used by this algorithm.
% the sintaxis "Bfs_fosmod(X,d,options)" is allowed in order to be
% similar to other Bfs_ functions.

if nargin==2;
    options = d;
end


m     = options.m;
show  = options.show;
[N,n] = size(X); % number of instances (samples), number of features
l     = zeros(m,1);
p     = ones(n,1);
q     = zeros(N,n);
for mm=1:m
    
    if mm==1
        Q = X;
    else
        Q = zeros(N,n);
        for j=1:n
            if p(j)
                alj = X(:,j);
                qm = alj;
                for k=1:mm-1
                    qm = qm - ((alj'*q(:,k))/(q(:,k)'*q(:,k)))*q(:,k);
                end
                Q(:,j) = qm;
            end
        end
        
    end
    
    
    C = zeros(n,n);
    for j=1:n
        if p(j)
            for i=1:n
                C(i,j) = Bfa_sqcorrcoef(X(:,i),Q(:,j));
            end
        end
    end
    Cm = mean(C);
    
    [i,l(mm)] = max(Cm);
    p(l(mm)) = 0;
    q(:,mm) = Q(:,l(mm));
    
    err = zeros(n,mm);
    for k=1:mm
        for j=1:n
            err(j,k) = Bfa_sqcorrcoef(X(:,j),q(:,k))*100;
        end
    end
    if mm>1
        serr = sum(err,2);
    else
        serr = err;
    end
    mserr = mean(serr);
    if show
        fprintf('%2d) selected feature=%4d error=%7.2f%%\n',mm,l(mm),100-mserr)
    end
end
selec = l;

