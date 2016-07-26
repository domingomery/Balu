% H = Bmv_homographySVD(m1,m2)
%
% Toolbox: Balu
%
%    Estimation of Homography Matrix using SVD decomposition.
%
%    m1 and m2 are n corresponding points in two views (m1(:,k) and m2(:,k)
%    are the k-th corresponding points (for k=1..n) stored as homogeneous
%    coordinates.
%
%    Example:
%       m1 = [rand(2,20);ones(1,20)];           % projection points in view 1 
%       H = rand(3,3);                          % original homography matrix
%       H = H/H(3,3)                            % normalization
%       m2 = H*m1;m2 = m2./(ones(3,1)*m2(3,:)); % projection points in view 2
%       Hs = Bmv_homographySVD(m1,m2);          % estimations of H
%       Hs = Hs/Hs(3,3)                         % normalization
%       m2s = Hs*m1;m2s = m2s./(ones(3,1)*m2s(3,:));
%       e = m2s-m2                              % reprojection error
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function H = Bmv_homographySVD(m1,m2)

np = size(m1,2);
xa = m1(1,:)';
ya = m1(2,:)';
xb = m2(1,:)';
yb = m2(2,:)';

A = zeros(2*np,9);

for i=1:np
    j = i*2-1;
    A(j:j+1,:) = [xa(i) ya(i) 1 0 0 0 -xa(i)*xb(i) -ya(i)*xb(i) -xb(i); 0 0 0 xa(i) ya(i) 1 -xa(i)*yb(i) -ya(i)*yb(i) -yb(i)];
end

[U,S,V] = svd(A); %#ok<ASGLU>
h = V(:,end);
H =[h(1) h(2) h(3); h(4) h(5)  h(6);h(7) h(8) h(9)];


