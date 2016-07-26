% [Xnew,a,b] = Bft_norm(X,normtype)
%
% Toolbox: Balu
%
%    Normalization of features X.
%    normtype = 1 for variance = 1 and mean = 0
%    normtype = 0 for max = 1, min = 0
%
%    Xnew = a*X + b
%
% Example:
%    load datareal
%    X = f(:,1:2);
%    figure(1); Bio_plotfeatures(X,d);
%    Xnew = Bft_norm(X,0); 
%    figure(2); Bio_plotfeatures(Xnew,d);
%
% (c) Grima, PUC-DCC, 2010: D. Mery, E.Cortazar
% http://dmery.ing.puc.cl

function [Xnew,a,b] = Bft_norm(X,normtype)

[N,M] = size(X);
if (normtype)
    mf    = mean(X);
    sf    = std(X);
    a     = ones(1,M)./sf;
    b     = -mf./sf;
else
    mi    = min(X);
    ma    = max(X);
    md    = ma-mi + (ma==mi);    
    a     = ones(1,M)./md;
    b     = -mi./md;
end
Xnew = X.*(ones(N,1)*a) + ones(N,1)*b;