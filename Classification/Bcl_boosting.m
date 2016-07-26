% ds      = Bcl_boosting(X,d,Xt,options)  Training & Testing together
% options = Bcl_boosting(X,d,options)     Training only
% ds      = Bcl_boosting(Xt,options)      Testing only
%
% Toolbox: Balu
%    Boosting classifier.
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.s portion of samples to be selected to design classifier C1.
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'boosting').
%       options.op1,op2,op3 contains training information about the three
%       classifiers used as weak classifiers.
%       options.dmin is min(d).
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%
%    See reference (Fig. 7):
%    Polikar, R. (2006): Ensemble Based Systems in Decision Making. IEEE
%    Circuits ans Systems Magazine, Third Quarter, 21-45.
%
%    Example: Training & Test together:
%       load datagauss                % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)         % plot feature space
%       op.s = 0.3;
%       ds = Bcl_boosting(X,d,Xt,op); % Boosting with s = 0.30
%       p = Bev_performance(ds,dt)    % performance on test data
%
%    Example: Training only
%       load datagauss                % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)         % plot feature space
%       op.s = 0.3;
%       op = Bcl_boosting(X,d,op);    % Boosting with s = 0.30
%
%    Example: Testing only (after training only example):
%       ds = Bcl_boosting(Xt,op);     % Boosting with s = 0.30
%       p = Bev_performance(ds,dt)    % performance on test data
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl
%

function [ds,options] = Bcl_boosting(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'boosting';
if train
    s = options.s;
    if s>1
        error('Bcl_boosting: s must be less than 1');
    end
    N = size(X,1);
    dmin = min(d);
    N1 = fix(s*N);
    d = d-dmin+1;
    op.p = [];

    % 1. Select N1 < N samples without replacement from X to create X1
    rn    = rand(N1,1);
    [i,j] = sort(rn);
    X1    = X(j,:);
    d1    = d(j);

    % 2. Tarining of C1 with X1
    op1 = Bcl_lda(X1,d1,op);
    ds1 = Bcl_lda(X,op1);

    % 3. Definition of dataset X2
    rn    = rand(N,1);
    [i,j] = sort(rn);
    Xr    = X(j,:);
    dr    = d(j);
    drs1  = ds1(j);

    i0 = find(drs1~=d);
    i1 = find(drs1==d);

    n2 = min([length(i0) length(i1)]);

    ii = [i0(1:n2);i1(1:n2)];

    X2 = Xr(ii,:);
    d2 = dr(ii);

    % 4. Training of classifier C2 using X2
    op2 = Bcl_lda(X2,d2,op);
    ds2 = Bcl_lda(X,op2);

    % 5. Definition of dataset X3
    ii = find(ds1~=ds2);
    if not(isempty(ii))
        X3 = Xr(ii,:);
        d3 = dr(ii);
        op3 = Bcl_lda(X3,d3,op);
    else
        op3 = [];
    end
    options.op1 = op1;
    options.op2 = op2;
    options.op3 = op3;
    options.dmin = dmin;
    ds = options;
end
if test
    op1  = options.op1;
    op2  = options.op2;
    op3  = options.op3;
    dmin = options.dmin;
    dt1  = Bcl_lda(Xt,op1);
    dt2  = Bcl_lda(Xt,op2);
    if not(isempty(op3))
        dt3 = Bcl_lda(Xt,op3);
        ds = dt1;
        ii = find(dt1~=dt2);
        if not(isempty(ii))
            ds(ii) = dt3(ii);
        end
%        ds = ds+dmin-1;
%    else
%        ds = dt1+dmin-1;
    end
    ds   = ds+dmin-1;
end

