% M = Bmv_reco3d2(m1,m2,A,B)
%
% Toolbox: Balu
%
% 3D reconstruction from n corresponding points
%
%    It returns a 3D point M that fullfils
%    the following projective equations:
% 
%    lambda1*m1 = A*M
%    lambda2*m2 = B*M
%    :
%    where mk are the 2D projection points of 3D point M
%    in image k; A and B are the corresponding 3x4 porjection
%    and lambdak are scale factors.
%
%    mk and M are given in homogeneous coordinates, i.e., mk
%    are 3x1 vectors, and M is a 4x1 vector. 
%
%
%    R. Hartley. A linear method for reconstruction from lines and
%    points. In 5th International Conference on Computer Vision 
%    (ICCV-95), pages 882-887, Cambridge, MA,1995.%
%
%    Example:
%       M = [1 2 3 1]'              % 3D point (X=1,Y=2,Z=3)
%       A = rand(3,4);              % proyection matrix for view 1
%       B = rand(3,4);              % proyection matrix for view 2
%       m1 = A*M; m1=m1/m1(3);      % proyection point in view 1
%       m2 = B*M; m2=m2/m2(3);      % proyection point in view 2
%       Ms = Bmv_reco3d2(m1,m2,A,B) % 3D reconstruction
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl
%

function M = Bmv_reco3d2(m1,m2,A,B)

H = [A;rand(1,4)];
while abs(det(H)) > 0.0001
   H = [A;rand(1,4)];
end
H1 = inv(H);

Bs = B*H1;
m2 = m2/m2(3);
m1 = m1/m1(3);
x2 = m2(1);
y2 = m2(2);
M = H1*[(y2*Bs(1,4)-x2*Bs(2,4))*m1/((x2*Bs(2,1:3)-y2*Bs(1,1:3))*m1); 1];
M = M/M(4);
