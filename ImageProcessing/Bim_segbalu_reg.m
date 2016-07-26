% [R,E,J] = Bim_segbalu(I,p)
%
% Toolbox: Balu
%  Segmentation of an object with homogeneous background.
% 
%  I: input image
%  p: threshold (default: p=-0.05) with p between -1 and 1. 
%      A positive value is used to dilate the segmentation, 
%      the negative to erode.
%  R: binary image of the object
%  E: binary image of the edge of the object
%  J: high contrast image of I.
% 
%  See details in:
%  Mery, D.; Pedreschi, F. (2005): Segmentation of Colour Food Images using 
%  a Robust Algorithm. Journal of Food Engineering 66(3): 353-360.
% 
%  Example:
%     I = imread('testimg1.jpg');
%     R = Bim_segbalu(I);
%     figure(1)
%     imshow(I); title('test image')
%     figure(2)
%     imshow(R); title('segmented image')
%
%     Repeat this examples for images testimg2, testimg3 and testimg4. Last
%     test image requires R = Bim_segbalu(I,-0.1) for better results.
% 
% See also Bim_segmowgli, Bim_segotsu, Bim_segkmeans, Bio_segshow.
%
% D.Mery, PUC-DCC, Apr. 2008-2010
% http://dmery.ing.puc.cl

function [R,E,J] = Bim_segbalu_reg(I,R,p)
ii = find(R==1);
Ir = I(:,:,1);
Ig = I(:,:,2);
Ib = I(:,:,3);
II = zeros(length(ii),1,3);
II(:,:,1) = Ir(ii);
II(:,:,2) = Ig(ii);
II(:,:,3) = Ib(ii);

Jo = Bim_rgb2hcm(double(II)/256);
J = zeros(size(R));
J(ii) = Jo;
t = graythresh(Jo);
if (~exist('p','var'))
    p = -0.05;
end
[R,E] = Bim_morphoreg_reg(J,R,t+p);