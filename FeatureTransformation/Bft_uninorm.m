% Xnew = Bft_uninorm(X)
%
% Toolbox: Balu
%
%    Normalization of features X: each row of Xnew has norm = 1
%
% Example:
%    load datareal
%    Xnew = Bft_uninorm(f);
%
% (c) Grima, PUC-DCC, 2013: D. Mery
% http://dmery.ing.puc.cl

function Xnew = Bft_uninorm(X)

[N,M] = size(X);
Xnew  = zeros(N,M);
for i=1:N
    Xnew(i,:) = X(i,:)/norm(X(i,:));
end
