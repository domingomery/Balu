% [M,err,ms] = Bmv_reco3dna(m,P)
%
% Toolbox: Balu
%
% 3D affine reconstruction from n corresponding points
%
%    It returns a 3D point M that fullfils
%    the following projective equations:
% 
%    m1 = P1*M
%    m2 = P2*M
%    :
%    where mk = m(:,k) are the 2D projection points of 3D point M
%    in image k; Pk = P(k*3-2:k*3,:) the corresponding 3x4 projection. The 
%    last row of a affine projection matrix must be [0 0 0 1]. factors
%    lambda k are 1 in affine projections.
%
%    mk and M are given in homogeneous coordinates, i.e., mk
%    are 3x1 vectors, and M is a 4x1 vector. 
%
%    ms is the reprojection of M in each view, ideally ms=m.
%    err is the reprojection error calculated as norm of ms-m.
%    The method used in this program is proposed in:
%
%    R. Hartley. A linear method for reconstruction from lines and
%    points. In 5th International Conference on Computer Vision 
%    (ICCV-95), pages 882-887, Cambridge, MA,1995.
%
%    Example:
%       M = [1 2 3 1]';                % 3D point (X=1,Y=2,Z=3)
%       P1 = [rand(2,4);0 0 0 1];      % proyection matrix for view 1
%       P2 = [rand(2,4);0 0 0 1];      % proyection matrix for view 2
%       P3 = [rand(2,4);0 0 0 1];      % proyection matrix for view 3
%       P4 = [rand(2,4);0 0 0 1];      % proyection matrix for view 4
%       m1 = P1*M;                     % proyection point in view 1
%       m2 = P2*M;                     % proyection point in view 2
%       m3 = P3*M;                     % proyection point in view 3
%       m4 = P4*M;                     % proyection point in view 4
%       P = [P1;P2;P3;P4];             % all projection matrices
%       m = [m1 m2 m3 m4];             % all proyection points
%       [Ms,err,ms] = Bmv_reco3dn(m,P) % 3D reconstruction
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [M,err,ms] = Bmv_reco3dna(m,P)

n = size(m,2); % number of images

Q = zeros(2*n,3);
r = zeros(2*n,1);

for k = 1:n
    p = P(k*3-2:k*3,:);
    Q(k*2-1:k*2,:) = -p(1:2,1:3); 
    r(k*2-1:k*2,:) = p(1:2,4) - m(1:2,k);
end

M = [(Q'*Q)\Q'*r; 1];

ms = zeros(3,n);
ms(:) = P*M;
d  = ms(1:2,:)-m(1:2,:);
err = sqrt(sum(d.*d,1));

