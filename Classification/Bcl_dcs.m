% ds      = Bcl_dcs(X,d,Xt,options)  Training & Testing together
% options = Bcl_dcs(X,d,options)     Training only
% ds      = Bcl_dcs(Xt,options)      Testing only
%
%   Toolbox: Balu
%      Dynamic classifier selection based on multiple
%             classifier behaviour after Giacinto (2001).
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%
%       obtions.b is a Balu classifier (i=1,...,n)
%          b(i).name    = classifier's name
%          b(i).options = options of classifier i
%
%       options.k   is the number of neighbours (e.g., k=7)
%       options.ths is the similarity-threshold (e.g., ths=0.8) between 0 and 1
%       options.thc is the classifier threshold (e.g., thc=0.7) between 0 and 1
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%
%
%    Giacinto, G. & Roli, F. (2001): Dynamic classifier selection based on
%    multiple classifier behaviour Pattern Recognition, 34: 1879 - 1881
%
%    Example Training & Testing
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
%       op.k   = 7;
%       op.ths = 0.8;
%       op.thc = 0.7;
%       ds = Bcl_dcs(X,d,Xt,op);             % DCS
%       p  = Bev_performance(ds,dt)          % performance on test data
%
%    Example: Training only
%       op = Bcl_dcs(X,d,op);                % DCS
%     
%    Example: Testing only (after training only example):
%       ds = Bcl_dcs(Xt,op);                 % DCS
%       p  = Bev_performance(ds,dt)          % performance on test data
%     
%    See also Bcl_ensemble.
%
%   D.Mery, PUC-DCC, 2010
%   http://dmery.ing.puc.cl

function [ds,options] = Bcl_dcs(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'dcs     ';
Bsclass = options.b;
op.b         = Bsclass;
op.tensemble = 'all';
op.show      = 0;
op.bs        = [];
kp  = options.k;
ths = options.ths;
thc = options.thc;
L = length(Bsclass);
if train
    op1 = Bcl_ensemble(X,d,op);
    options.X = X;
    options.d = d;
    options.op1 = op1;
    ds = options;
end
if test
    X = options.X;
    d = options.d;
    op1 = options.op1;
    nt  = size(Xt,1);
    ds = zeros(nt,1);
    n1  = size(X,1);
    DS1 = Bcl_ensemble(X ,op1);
    DS2 = Bcl_ensemble(Xt,op1);
    for j=1:nt
        dX    = X-ones(n1,1)*Xt(j,:);
        dd2   = sum(dX.*dX,2);
        [I,J] = sort(dd2);
        K     = J(1:kp);
        Ds    = d(K);
        D1s   = DS1(K,:);
        D2s   = DS2(j,:);
        S     = mean(D1s == ones(kp,1)*D2s,2);
        ii    = find(S>ths);
        if not(isempty(ii))
            CLA   = sum(D1s(ii,:) == (Ds(ii)*ones(1,L)),1)/length(ii);
            [I,J] = max(CLA);
            if I>=thc
                ds(j) = D2s(J);
            else
                ds(j) = round(mean(D2s));
            end
        else
            ds(j) = round(mean(D2s));
        end
    end
end

