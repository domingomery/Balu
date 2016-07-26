% function m = Bmv_lines2point(l1,l2)
%
% Toolbox Balu
%
% 2D point m computed as the intersection of two 2D lines (l1 and l2).
%
% m = Bmv_lines2point(l1,l2) returns the 2D point m defined
%    as the intersection of lines l1 and l2.
%    m, l1, l2 are 3x1 homogeneous vectors. The result
%    m is given as m = [x y 1]' where m'*1l = m'*l2 = 0
%
% Example:
%    l1 = [2 1 1]; l2 = [1 2 1];
%    m = Bmv_lines2point(l1,l2)
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function m = Bmv_lines2point(l1,l2)
m = cross(l1,l2);
m = m/m(3);