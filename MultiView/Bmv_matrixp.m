% function P = Bmv_matrixp(f)
%
% Toolbox Balu
%
% Perspective proyection matrix 3D->2D.
%
% It returns the 3x4 perspective proyection matrix
%    depending on focal distance f.
%
%    Bmv_matrixp(f) is equal to 
%    [f 0 0 0
%     0 f 0 0
%     0 0 1 0]
%
% Example:
%
%     f = 10;              % focal distance
%     P = Bmv_matrixp(f);  % Perspective matrix
%     M = [1 2 3 1]';      % 3D point
%     m = P*M; m = m/m(3)  % 2D projected point
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl
%
function P = Bmv_matrixp(f)
P = [f 0 0 0;0 f 0 0;0 0 1 0];
