% [T,U,P,Q,W,B] = Bft_plsr(X,d,options)
%
% Toolbox: Balu
%    Feature transformation using Partial Least Squares Regression with
%    NIPALS algorithm.
%    X: Input matrix with features
%    d: Vector with ideal classifcation.
%    m: Number of principal components to be selected.
%    T: Loadings of X (m transformed features). Matrix T is Xo*W', where
%       Xo is normalized X and W corresponds to the weights.
%    U: Loadings of Y (m transformed features). Matrix U is T*B.
%    P: Principal components of X
%    Q: Principal components of Y (defined from d, one binary column
%       per class)
%    W: X weights
%    B: A vector of regression coefficients
%
%    This code was adapted from nipals.m developed by
%    http://www.cdpcenter.org/files/plsr after Geladi (1986)
%
%    Reference:
%    Geladi, P. & Kowalski, B. R. (1986): Partial least-squares regression: a
%    tutorial Analytica Chimica Acta, (185):1-17.
%
%    Example:
%       load datareal
%       s = [279 235 268 230 175 165 207 160 269 157]; %indices using Example
%                                                      %of Bfs_sfs.
%       X = f(:,s);
%       T = Bft_plsr(X,d,6);
%       op.p = [];
%       ds1 = Bcl_lda(T,d,T,op);
%       Bev_performance(d,ds1)
%
%       % Comparison with first 6 SFS features
%       ds2 = Bcl_lda(X(:,1:6),d,X(:,1:6),op);
%       Bev_performance(d,ds2)
%
% D.Mery, PUC-DCC, May. 2010
% http://dmery.ing.puc.cl
%

function [T,U,P,Q,W,B] = Bft_plsr(X,d,options)

if isfield(options,'m')
    m = options.m;
    if isfield(options,'show')
        show = options.show;
    end
else
    m = options;
    show = 0;
end


a = m;
[n,m] = size(X);           % The number of observations and inputs


Y = d2Y(d);
p = size(Y,2);             % The number of outputs

P = zeros(m, a, 'single'); % The principal components of _X_
Q = zeros(p, a, 'single'); % The principal components of _Y_
W = zeros(m, a, 'single'); % The _X_ weights
B = zeros(a, a, 'single'); % A vector of regression coefficients
T = zeros(n, a, 'single'); % The loadings of _X_
U = zeros(n, a, 'single'); % The loadings of _Y_

X = Bft_norm(X,1);         % normalization: mean(X) = 0, std(X) = 1
Y = Bft_norm(Y,1);         %                mean(Y) = 0, std(Y) = 1

Ye = Y;

if show
    ff = Bio_statusbar('PLSR progress');
end

for h=1:a
    if show
        ff = Bio_statusbar(h/a,ff);
    end
    t = zeros(n, 1, 'single');
    u = Ye(:,1);
    err = 1;
    while err>1e-4
        w   = X'*u;
        w   = w/norm(w);
        tn  = X*w;
        q   = Ye'*tn;
        q   = q/norm(q);
        u   = Ye*q;
        err = norm(t - tn)/norm(tn);
        t   = tn;
    end
    
    p = X'*t;
    p = p/norm(p);
    P(:,h) = p;
    Q(:,h) = q;
    W(:,h) = w;
    T(:,h) = t;
    U(:,h) = u;
    % One difference from Geladi's algorithm is that
    % we never calculate the X-residual, so we need
    % to calculate the regression coefficients a
    % bit differently. (See Eq. 14 in Geladi.)
    B(1:h,h) = (T(:,1:h)'*T(:,1:h)) \ T(:,1:h)'*u;
    Ys = T*B*Q';
    Ye = Y - Ys;
end
if show
    delete(ff)
end