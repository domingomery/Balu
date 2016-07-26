% J = Bim_lin(I,t)
%
% Toolbox: Balu
%    Lineal enhancement of image I from 0 to 255.
%
%    Input data:
%       I grayvalue image.
%
%    Output:
%       J: enhanced image so that
%       J = m*I + b; where min(J) = 0, max(J) = 255.
%       J is uint8.
%
%    Example:
%       X = imread('testimg2.jpg');
%       I = rgb2gray(X);
%       figure(1)
%       imshow(I); title('original image')
%       J = Bim_lin(I);
%       figure(2)
%       imshow(J); title('transformed image')
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function J = Bim_lin(I,t)

if ~exist('t','var')
    t = 0;
end

I = double(I);
mi = min(I(:));
mx = max(I(:));
d  = mx - mi;
if d==0
    J = I/mx*128;
else


    if t==0
        J = (I-mi)/d*255;
    else
        mx = max(abs([mi mx]));
        J  = I*127/mx+128;        
    end

end
J = uint8(round(J));