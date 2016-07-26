% F = Bmv_trifocalSVD(m1,m2,m3)
%
% Toolbox: Balu
%
%    Estimation of Trifocal Tensors using SVD decomposition.
%
%    m1, m2 and m3 are n corresponding points in 3 views m1(:,k), m2(:,k)
%    and m3(:,k) are the k-th corresponding points (for k=1..n) stored as 
%    homogeneous coordinates.
%
%    Example:
%       A = rand(3,4);                         % Projection matrix for view 1
%       B = rand(3,4);                         % Projection matrix for view 2
%       C = rand(3,4);                         % Projection matrix for view 3 
%       n = 20;                                % 20 corresponding points
%       M = [rand(3,n);ones(1,n)];             % n 3D points
%       m1 = A*M;m1 = m1./(ones(3,1)*m1(3,:)); % projection points in view 1
%       m2 = B*M;m2 = m2./(ones(3,1)*m2(3,:)); % projection points in view 2
%       m3 = C*M;m3 = m3./(ones(3,1)*m3(3,:)); % projection points in view 3
%       T = Bmv_trifocalSVD(m1,m2,m3);         % Trifocal tensors
%       m3s = Bmv_reproj3(m1,m2,T);            % Reprojection of points in view 3
%       [m3s-m3]                               % Difference must be zero
%       
%
% (c) D.Mery, PUC-DCC, 2012
% http://dmery.ing.puc.cl


function T = Bmv_trifocalSVD(m1,m2,m3)
n1 = size(m1,2);
n2 = size(m2,2);
n3 = size(m3,2);
if (n1==n2) && (n2==n3)
    n = n1;
else
    error('m1, m2 and m3 have different number of points');
end

A = zeros(4*n,27);

for i = 1 : n
    
    q   = 4*(i-1);
    m1i = m1(:,i);
    x2i = m2(1,i );
    y2i = m2(2,i );
    x3i = m3(1,i );
    y3i = m3(2,i );
        
    % Trilinearity D1
    A(q+1,1:3)   = -m1i;
    A(q+1,7:9)   =  m1i*x3i;
    A(q+1,19:21) =  m1i*x2i;
    A(q+1,25:27) = -m1i*x2i*x3i;

    % Trilinearity D2
    A(q+2,4:6)   = -m1i;
    A(q+2,7:9)   =  m1i*y3i;
    A(q+2,22:24) =  m1i*x2i;
    A(q+2,25:27) = -m1i*x2i*y3i;

    % Trilinearity D3
    A(q+3,10:12) = -m1i;
    A(q+3,16:18) =  m1i*x3i;
    A(q+3,19:21) =  m1i*y2i;
    A(q+3,25:27) = -m1i*x3i*y2i;

    % Trilinearity D4
    A(q+4,13:15) = -m1i;
    A(q+4,16:18) =  m1i*y3i;
    A(q+4,22:24) =  m1i*y2i;
    A(q+4,25:27) = -m1i*y2i*y3i;
    
end

% ||Ax|| -> min s.t. ||x||=1

[U,S,V] = svd(A);

x = V(:,end);

T = zeros(3,3,3);

% T(1,1,1) = x(1);
% T(2,1,1) = x(2);
% T(3,1,1) = x(3);
% T(1,2,1) = x(10);
% T(2,2,1) = x(11);
% T(3,2,1) = x(12);
% T(1,3,1) = x(19);
% T(2,3,1) = x(20);
% T(3,3,1) = x(21);
% 
% T(1,1,2) = x(4);
% T(2,1,2) = x(5);
% T(3,1,2) = x(6);
% T(1,2,2) = x(13);
% T(2,2,2) = x(14);
% T(3,2,2) = x(15);
% T(1,3,2) = x(22);
% T(2,3,2) = x(23);
% T(3,3,2) = x(24);
% 
% T(1,1,3) = x(7);
% T(2,1,3) = x(8);
% T(3,1,3) = x(9);
% T(1,2,3) = x(16);
% T(2,2,3) = x(17);
% T(3,2,3) = x(18);
% T(1,3,3) = x(25);
% T(2,3,3) = x(26);
% T(3,3,3) = x(27);

T(:) = x([1 2 3 10 11 12 19 20 21 4 5 6 13 14 15 22 23 24 7 8 9 16 17 18 25 26 27]);
