% function m3s = Bmv_reproj3(m1,m2,T,method)
%
% Toolbox Balu:
%
% Reprojection of point m3 from m1, m2 and trifocal tensors
%
% m3s = reproj3(m1,m2,T) returns the reprojection of m3 from
%    corresponding points m1 and m2 in image 1 and 2 respectivelly
%    using trifocal tensors.
%    method = 1 uses the first two trilinearities
%    method = 2 uses all four trilinearities (default)
%
% Example:
%
%       A   = rand(3,4);           % Projection matrix for view 1
%       B   = rand(3,4);           % Projection matrix for view 2
%       C   = rand(3,4);           % Projection matrix for view 3
%       T   = Bmv_trifocal(A,B,C); % Trifocal tensors
%       M   = rand(4,1);           % 3D points
%       m1  = A*M;                 % 2D point in view 1 (third term can be ~=1)
%       m2  = B*M;                 % 2D point in view 2 (third term can be ~=1)
%       m3s = Bmv_reproj3(m1,m2,T);
%       m3  = C*M; m3 = m3/m3(3);  % 2D point in view 3 
%       [m3 m3s]
%
% (c) D.Mery, PUC-DCC, 2012
% http://dmery.ing.puc.cl
%
function m3s = Bmv_reproj3(m1,m2,T,method)

m1 = m1./(ones(3,1)*m1(3,:));
m2 = m2./(ones(3,1)*m2(3,:));



if ~exist('method','var')
    method = 2;
end

n = size(m1,2);

m3s = zeros(3,n);

for k=1:n
    m3s(:,k) =  [T(:,1,1)-m2(1,k)*T(:,3,1)  T(:,1,2)-m2(1,k)*T(:,3,2) T(:,1,3)-m2(1,k)*T(:,3,3)]'*m1(:,k);
end

m3s = m3s./(ones(3,1)*m3s(3,:));

if method==2
    m3s2 = zeros(3,n);
    for k=1:n
        m3s2(:,k) =  [T(:,2,1)-m2(2,k)*T(:,3,1)  T(:,2,2)-m2(2,k)*T(:,3,2) T(:,2,3)-m2(2,k)*T(:,3,3)]'*m1(:,k);
    end
    m3s2 = m3s2./(ones(3,1)*m3s2(3,:));
    m3s = (m3s+m3s2)/2;
end

