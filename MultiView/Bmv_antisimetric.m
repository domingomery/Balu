% U = Bmv_antisimetric(u)
%
% Toolbox: Balu
%
%    Antisimetric matrix
%
%    antisimetric(u) returns the antisimetric matrix of a
%    a 3x1 vector u.
%
%    U = antisimetric(u) is a 3x3 matrix that
%    U*v = cross(u,v) where cross(u,v) is the cross
%    product between u and a 3x1 vector v.
%
%    Example:
%       u = [1 2 3]'
%       v = rand(3,1)
%       w1 = cross(u,v)
%       U  = Bmv_antisimetric(u)
%       w2 = U*v              % w1 must be equal to w2
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function U = Bmv_antisimetric(u)

U = [  0   -u(3)  u(2)
      u(3)   0   -u(1)
     -u(2)  u(1)   0  ];