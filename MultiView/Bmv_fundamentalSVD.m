% F = Bmv_fundamentalSVD(m1,m2)
%
% Toolbox: Balu
%
%    Estimation of Fundamental Matrix using SVD decomposition (Hartley, 
%    p.281) with det(F) = 0, norm(F(:)) = 1.
%
%    m1 and m2 are n corresponding points in two views (m1(:,k) and m2(:,k) 
%    are the k-th corresponding points (for k=1..n) stored as homogeneous
%    coordinates.
%
%    Example:
%       A = rand(3,4);                         % Projection matrix for view 1
%       B = rand(3,4);                         % Projection matrix for view 2
%       n = 20;                                % 20 corresponding points
%       M = [rand(3,n);ones(1,n)];             % n 3D points
%       m1 = A*M;m1 = m1./(ones(3,1)*m1(3,:)); % projection points in view 1
%       m2 = B*M;m2 = m2./(ones(3,1)*m2(3,:)); % projection points in view 2
%       F = Bmv_fundamentalSVD(m1,m2);         % Fundamental matrix using tensors
%       diag(m2'*F*m1)                         % epipolar constraint must be zero
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function F = Bmv_fundamentalSVD(m1,m2)

n = size(m1,2);
x1 = m1(1,:)';
y1 = m1(2,:)';
x2 = m2(1,:)';
y2 = m2(2,:)';


A = [x2.*x1 x2.*y1 x2 y2.*x1 y2.*y1 y2 x1 y1 ones(n,1)];


[U,S,V] = svd(A);

f = V(:,end);

Fs = zeros(3,3);
Fs(:) = f;
Fs = Fs';

[U,S,V] = svd(Fs);
S(3,3) = 0;
F = U*S*V';