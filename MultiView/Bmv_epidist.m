% d = Bmv_epidist(m1,m2,F,method)
%
% Toolbox: Balu
%
% Distance from m2 to epipolar line l2 = F*m1
%
%    d = distance2(m1,m2,F,'method') returns the distance 
%    error. The posible corresponding points are m2 and m1.
%    F is the fundamental matrix. The distance is calculated 
%    using the following methods:
%   
%    method = 'euclidean': uses Euclidean distance from m2 to l=F*m1.
%
%    method = 'sampson': uses Sampson distance.
%
%    If no method is given, 'euclidean' will be assumed as default.
%
%    Both methods can be found in:
%
%    R. Hartley and A. Zisserman. Multiple View Geometry in Computer
%    Vision. Cambridge University Press, 2000.
%
%    Example:
%       A = rand(3,4);                      % Projection matrix for view 1
%       B = rand(3,4);                      % Projection matrix for view 1
%       M = [1 2 3 1]';                     % 3D point (X=1,Y=2,Z=3)
%       m1 = A*M;m1 = m1/m1(3);             % projection point in view 1
%       m2 = B*M;m2 = m2/m2(3);             % projection point in view 2
%       F = Bmv_fundamental(A,B,'tensor');  % Fundamental matrix using tensors
%       d = Bmv_epidist(m1,m2,F,'sampson')  % sampson distance to epipolar line
%
%    See also Bmv_epiplot.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function d = Bmv_epidist(m1,m2,F,method)

if ~exist('method','var')
   method = 'euclidean';
end

d0 = abs(m2'*F*m1);

switch lower(method)
   case 'sampson'
      l1 = F*m1;
      l2 = F'*m2;
      d = d0/sqrt(l1(1)^2+l1(2)^2+l2(1)^2+l2(2)^2);
    otherwise % euclidean
      l = F*m1;
      d = d0/sqrt(l(1)^2+l(2)^2);
end
