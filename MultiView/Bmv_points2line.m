% function l = Bmv_points2line(m1,m2)
%
% Toolbox Balu
%
% 2D line l that contains two 2D points (m1 and m2).
%
% l = points2line(m1,m2) returns the 2D line l defined as the line that 
%     contains the 2D points m1 and m2. l, m1, m2 are 3x1 homogeneous 
%     vectors. Points m1 and m2 can be 2x1 vectors. The result l is given 
%     as  [a b c]' where m1'*l = m2'*l = 0, i.e., the equation of the line
%     is a*x + b*y + c = 0
%
% Example:
%    m1 = [2 1]; m2 = [1 2];
%    l = Bmv_points2line(m1,m2)
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function l = Bmv_points2line(m1,m2)
if length(m1)==2
    mm1 = [m1(1) m1(2) 1]';
else
    mm1 = m1;
end
if length(m2)==2
    mm2 = [m2(1) m2(2) 1]';
else
    mm2 = m2;
end
l = cross(mm1,mm2);
