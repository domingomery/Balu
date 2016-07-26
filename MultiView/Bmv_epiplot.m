% ell = Bmv_epiplot(F,m1)
%
% Toolbox: Balu
%
%    Plot of epipolar line.
% 
%    The epipolar line is ell = F*m1. F is the Fundamental Matrix and m1 is
%    a point image 1 in homogeneous coordinates.
%
%    Example:
%       I1 = imread('testimg5.jpg');            % Image 1
%       figure(1);imshow(I1); hold on
%       I2 = imread('testimg6.jpg');            % Image 2
%       figure(2);imshow(I2); hold on
%       F  = Bmv_fundamentalSIFT(I1,I2);        % F matrix estimation
%       while(1)
%          figure(1);
%          disp('click a point in Figure 1...') % click
%          p = vl_click; m1 = [p(2) p(1) 1]';
%          plot(p(1),p(2),'g+')
%          figure(2)
%          Bmv_epiplot(F,m1)
%       end
%
%    See also Bmv_fundamentalSIFT, Bmv_fundamentalRANSAC, Bmv_fundamentalSVD. 
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function ell = Bmv_epiplot(F,m1,col)
if ~exist('col','var')
    col = 'b';
end
ell = F*m1;
ax = axis;
x = ax(3:4)';
a  = ell(1);
b  = ell(2);
c  = ell(3);
y = -(c+a*x)/b;
plot(y,x,col)
drawnow