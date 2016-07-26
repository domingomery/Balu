% [F,inliers] = Bmv_fundamentalRANSAC(m1,m2)
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
%       A = rand(3,4);                          % Projection matrix for view 1
%       B = rand(3,4);                          % Projection matrix for view 2
%       n = 20;                                 % 20 corresponding points
%       M = [rand(3,n);ones(1,n)];              % n 3D points
%       m1 = [A*M rand(3,6)]; m1 = m1./(ones(3,1)*m1(3,:));  % projection points in view 1 plus 6 outliers 
%       m2 = [B*M rand(3,6)] ;m2 = m2./(ones(3,1)*m2(3,:));  % projection points in view 2 plus 6 outliers 
%       options.dmax = 0.1;
%       options.iter = 1000;
%       [F,i] = Bmv_fundamentalRANSAC(m1,m2);   % Fundamental matrix using tensors
%       diag(m2'*F*m1)                          % epipolar constraint near to zero
%                                               % for the first 20 points.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [F,inliers] = Bmv_fundamentalRANSAC(m1,m2,options)

if nargin==2
    options.iter = 1000;
    options.dmax = 3;
end
if ~isfield(options,'iter');
   options.iter = 1000;
end
if ~isfield(options,'dmax');
   options.dmax = 3;
end
iter = options.iter;
dmax = options.dmax;

[mm,n]   = size(m1);
[mm2,n2] = size(m2);
if n~=n2
    error('Bmv_fundamentalRANSAC: number of points in m1 and m2 must be equal.');
end

if mm==2
    m1 = [m1;ones(1,n)];
end

if mm2==2
    m2 = [m2;ones(1,n)];
end

if or(size(m2,1)~=3,size(m1,1)~=3)
    error('Bmv_fundamentalRANSAC: points m1 and m2 must be bidimensional.');
end

x1 = m1(1,:)';
y1 = m1(2,:)';
x2 = m2(1,:)';
y2 = m2(2,:)';

A = [x2.*x1 x2.*y1 x2 y2.*x1 y2.*y1 y2 x1 y1 ones(n,1)];

outliersmax = Inf;
for i=1:iter
    r = rand(n,1);
    [ri,rj] = sort(r); 
    rj = rj(1:10);
    m1r = m1(:,rj);
    m2r = m2(:,rj);
    Fr = Bmv_fundamentalSVD(m1r,m2r);
    lr = (Fr*m1)';
    lrs = lr.*lr;
    lrs = sqrt(lrs(:,1)+lrs(:,2));
    Frt = Fr';
    d = abs(A*Frt(:))./lrs;
    outliers = sum(d>dmax);
    if outliers<outliersmax
        outliersmax = outliers;
        F = Fr;
        inliers = d<=dmax;
    end
end
