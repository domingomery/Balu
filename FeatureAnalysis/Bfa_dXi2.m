% d = Bfa_dXi2(X,Y)
%
% Toolbox: Balu
%
%    Xi^2 distance between two vectors X and Y
%
%    Example:
%       X = (1:50)';
%       Y = X + randn(50,1);
%       d = Bfa_dXi2(X,Y)
%    
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function d = Bfa_dXi2(X,Y)
X     = double(X);
Y     = double(Y);
s     = X+Y;
d2    = (X-Y).*(X-Y);
d     = sum(d2(s>0)./s(s>0));
