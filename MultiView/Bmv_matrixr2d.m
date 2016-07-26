% function R = Bmv_matrixr2d(theta)
%
% Toolbox Balu:
%
% 2D rotation matrix.
%
% It returns the 2D rotation matrix given by a rotation
%    theta (given in radians):
%
%  R = [ cos(theta) -sin(theta)
%      sin(theta)  cos(theta)];
%
%
% Example:
%
%    R = matrixr2d(pi/3) 
%
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function R = Bmv_matrixr2d(theta)

R = [ cos(theta) -sin(theta)
      sin(theta)  cos(theta)];

