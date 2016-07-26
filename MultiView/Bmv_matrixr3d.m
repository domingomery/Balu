% function R = Bmv_matrixr3d(wx,wy,wz)
%
% Toolbox Balu:
%
% 3D rotation matrix.
%
% It returns the 3D rotation matrix given by a rotation
%    arround z, y and x axes where the rotation angles are wz, wy, and
%    wx respectively. The angles are given in radians.
%
%    R = Bmv_matrixr3d(wx,wy,wz) is equal to Rx*Ry*Rz where 
%    
%
%
%    Rz = [  cos(wz)  sin(wz)    0
%           -sin(wz)  cos(wz)    0
%             0         0       1];
%
%    Ry = [ cos(wy)     0     -sin(wy)
%              0        1       0
%           sin(wy)     0      cos(wy)]; 
%
%    Rx = [    1        0        0
%              0      cos(wx)  sin(wx)
%              0     -sin(wx)  cos(wx)];
%
%  Example:
%     R = Bmv_matrixr3d(pi/2,-pi,0)
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function R = Bmv_matrixr3d(wx,wy,wz)

R = [
cos(wy)*cos(wz)                          cos(wy)*sin(wz)                         -sin(wy) 
sin(wx)*sin(wy)*cos(wz)-cos(wx)*sin(wz)  sin(wx)*sin(wy)*sin(wz)+cos(wx)*cos(wz)  sin(wx)*cos(wy) 
cos(wx)*sin(wy)*cos(wz)+sin(wx)*sin(wz)  cos(wx)*sin(wy)*sin(wz)-sin(wx)*cos(wz)  cos(wx)*cos(wy)];


