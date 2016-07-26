% selec = Bfs_rank(X,d,options)%
%
% Toolbox: Balu
%    Feature selection based on command rankfeatures (from MATLAB 
%    Bioinformatics Toolbox) that ranks ranks key features by class
%    separability criteria.
%
%    input: X feature matrix
%           options.m number of features to be selected
%           options.criterion can be:
%      'ttest' (default) Absolute value two-sample T-test with pooled
%                        variance estimate
%      'entropy'         Relative entropy, also known as Kullback-Lieber
%                        distance or divergence
%      'brattacharyya'   Minimum attainable classification error or
%                        Chernoff bound
%      'roc'             Area between the empirical receiver operating
%                        characteristic (ROC) curve and the random classifier
%                        slope
%      'wilcoxon'        Absolute value of the u-statistic of a two-sample
%                        unpaired Wilcoxon test, also known as Mann-Whitney
%      
%     Notes: 1) 'ttest', 'entropy', and 'brattacharyya' assume normal
%     distributed classes while 'roc' and 'wilcoxon' are nonparametric tests,
%     2) all tests are feature independent.
%
%    output: selec selected features
%
% Example:
%    load datareal
%    op.m = 10;                     % 10 features will be selected
%    op.criterion = 'roc';          % ROC criterion will be used
%    op.show = 1;                   % display results
%    s = Bfs_rank(f,d,op);          % index of selected features
%    X = f(:,s);                    % selected features
%    Xn = fn(s,:)                   % list of feature names
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl
%

function selec = Bfs_rank(X,d,options)

m = options.m;
criterion = options.criterion;

dmin = min(d);
dmax = max(d);
k = dmax-dmin+1;

if not(exist('criterion','var'))
    criterion = 'ttest';
end

if k<2
    error('Bfs_rank: Number of classes of d must be greater than 1.');
end

if ~exist('rankfeatures','file')
    error('Bfs_rank: This function requires Bioinformatics Toolbox.');
end



if k==2;
    idx = rankfeatures(X',d,'criterion',criterion);
    selec = idx(1:m);
else
    [N,M] = size(X);
    S = zeros(M,k);
    for j=1:k
        dj = (d==j)+1;
        S(:,j) = rankfeatures(X',dj,'criterion',criterion);
    end    
    T = S';
    idx = T(:);
    i = 1;
    selec = zeros(m,1);
    selec(1) = idx(1);
    j=2;
    while i<=m
        if not(ismember(idx(j),selec))
            selec(i) = idx(j);
            i = i+1;
        end
        j = j + 1;
    end    
end


        
        
        
