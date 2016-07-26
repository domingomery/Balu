% function [T,P,B,Y] = Bft_pcr(X,d,m)
%
% Toolbox: Balu
%
%    Principal Component Regression. 
%
%    X: Input matrix with features
%    d: Vector with ideal classifcation.
%    m: Number of principal components to be selected.
%    T: Loadings of X (m transformed features). Matrix T is Xo*W', where
%       Xo is normalized X and W corresponds to the wights.
%    P: Principal components of X
%    B: A vector of regression coefficients
%    Y: A NxC matrix Y(i,j) = 1 if sample i belongs to class j
%       N is the number of sanples, C the number of classes.
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
%       [T,P,B,Y] = Bft_pcr(X,d,6);
%       Ys = T*B;
%       mse(Ys-Y)                                      % mean square error
%
%    
% D.Mery, PUC-DCC, May. 2010
% http://dmery.ing.puc.cl
%

function [T,P,B,Y] = Bft_pcr(X,d,m)

N = size(X,1);           % The number of observations and inputs

dmin = min(d);
dmax = max(d);
Y = zeros(N,dmax-dmin+1);
for i=dmin:dmax
    Y(:,i) = d==i;
end

X = Bft_norm(X,1);
[T,lambda,P] = Bft_pca(X,m);
B = (T'*T)\T'*Y;
