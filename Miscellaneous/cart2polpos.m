% [t,r] = cart2polpos(x,y,s)
% 
% Toolbox: Balu
%    transform Cartesian to polar coordinates. The angle t will be always
%    positive. If s==0, the t is between 0 and 2*pi. If s > 0, the angles are
%    subdivided into s bins, and t means in which bin the angle is located.
%
% Examples:
%  [t,r] = cart2polpos(1,1)  % it is the same as [t,r] = cart2polpos(1,1,0)
%  the result for t is pi/4 and for r = sqrt(2)/2
%
%  [t,r] = cart2polpos(1,-1)  % it is the same as [t,r] = cart2polpos(1,1,0)
%  the result for t is 7*pi/4 and for r = sqrt(2)/2
%
%  [t,r] = cart2polpos(1,-1,4) % example that computes the quadrant
%  the result for t is 4 and for r = sqrt(2)/2
%
% D.Mery, PUC-DCC, 2014
% http://dmery.ing.puc.cl

function [t,r] = cart2polpos(x,y,s)

if ~exist('s','var')
    s=0;
end

[t,r] = cart2pol(x,y);
t(t<0) = t(t<0)+2*pi;
if s>0
    t = fix(t/2/pi*s)+1;
end
