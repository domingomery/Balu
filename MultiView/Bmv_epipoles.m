% [e1,e2] = Bmv_epipoles(F)
%
% Toolbox: Balu
%
% Epipoles of a two-view system from fundamental matrix.
%
% [e1,e2] = epipoles(F) returns:
%    e1: epipole in view 1
%    e2: epipole in view 2
%    where e1 and e2 are 3x1 homogeneous vectors.
%    e1(3) = e2(3) is 1.
%
%    Example:
%       A = rand(3,4);                      % Projection matrix for view 1
%       B = rand(3,4);                      % Projection matrix for view 1
%       F = Bmv_fundamental(A,B,'pseudo');  % Fundamental matrix
%       [e1,e2] = Bmv_epipoles(F);          % Epipoles
%       m1 = rand(3,1);
%       m2 = rand(3,1)
%       m2'*F*e1                            % epipol e1 belongs to epipolar line l1
%       e2'*F*m1                            % epipol e1 belongs to epipolar line l2
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [e1,e2] = Bmv_epipoles(F)

e1 = null(F);
if isempty(e1)
    error('Error: determinant of F is not equal to zero (|F|=%f).')
else
    e1 = null(F);
    e1 = e1/e1(3);
    e2 = null(F');
    e2 = e2/e2(3);
end
