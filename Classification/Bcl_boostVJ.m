% ds      = Bcl_boostVJ(X,d,Xt,options)  Training & Testing together
% options = Bcl_boostVJ(X,d,options)     Training only
% ds      = Bcl_boostVJ(Xt,options)      Testing only
%
% Toolbox: Balu
%    Boosting algorithm after Viola Jones. It uses only one feature
%    per weak classifier.
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.iter is the number of iterations (default=25).
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.dmin contains min(d).
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'boostVJ ').
%       options.at,pt,tht,jt,sat are training parameters used by testing.
%       options.jt are the selected features.
%       options.tht are the thresholds for each features
%       options.pt(t) = 1 means feature(jt(t))<tht(t)
%
%    See reference:
%    Viola, P., Jones, M. (2004): Robust real-time ob ject detection.
%    International Journal of Computer Vision (57):137-154.
%
%    Example: Training & Test together:
%       load datareal                           % real data (200 samples, 279 features,
%                                               % 2 classes)
%       op.iter = 20;
%       [X1,d1,X2,d2] = Bds_stratify(f,d,0.75); % Training and test data
%       ds = Bcl_boostVJ(X1,d1,X2,op);          % Bviola with 20 iterarions
%       p = Bev_performance(ds,d2)              % performance on test data
%
%    Example: Training only
%       load datareal                           % real data (200 samples, 279 features,
%                                               % 2 classes)
%       op.iter = 20;
%       [X1,d1,X2,d2] = Bds_stratify(f,d,0.75); % Training and test data
%       op = Bcl_boostVJ(X1,d1,op);             % Bviola with 20 iterarions
%
%    Example: Testing only (after training only example):
%       ds = Bcl_boostVJ(X2,op);                % Bviola with 20 iterarions
%       p = Bev_performance(ds,d2)              % performance on test data
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


function [ds,options] = Bcl_boostVJ(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'boostVJ ';
if train
    T = options.iter;

    dmin = min(d);
    dmax = max(d);
    N = dmax-dmin+1;
    if N~=2
        error('Bcl_viola works only with two classes. Data has %d classes.',N);
    end

    d = d-dmin;
    [n,m] = size(X);
    i0 = find(d==0);
    i1 = find(d==1);
    n0 = length(i0);
    n1 = length(i1);

    w = zeros(n,1);
    w(i0) = 1/2/n0;
    w(i1) = 1/2/n1;

    sat = 0;
    at  = zeros(T,1);
    pt  = zeros(T,1);
    jt  = zeros(T,1);
    tht = zeros(T,1);
    for t = 1:T
        w = w/sum(w);
        et = Inf;
        for j = 1:m
            [hj,op] = Bcl_weakc(X(:,j),d,X(:,j),[]);
            pj  = op.p;
            thj = op.th;
            ej = w'*abs(hj-d);
            if ej<et
                et = ej;
                jj = j;
                pp = pj;
                tt = thj;
                ht = hj;
            end
        end
        bt     = et/(1-et);
        at(t)  = log(1/bt);
        pt(t)  = pp;
        jt(t)  = jj;
        tht(t) = tt;
        sat = sat + at(t);
        e   = abs(ht-d);
        w  = w.*bt.^(1-e);
    end
    options.at   = at;
    options.pt   = pt;
    options.tht  = tht;
    options.jt   = jt;
    options.dmin = dmin;
    options.sat  = sat;
    ds = options;
end
if test
    at   = options.at;
    pt   = options.pt;
    tht  = options.tht;
    jt   = options.jt;
    dmin = options.dmin;
    sat  = options.sat;
    nt   = size(Xt,1);
    sht  = zeros(nt,1);
    for t = 1:options.iter
        sht = sht+at(t)*(pt(t)*Xt(:,jt(t))<pt(t)*tht(t));
    end
    ds = sht>=0.5*sat;
    ds = ds+dmin;
end


