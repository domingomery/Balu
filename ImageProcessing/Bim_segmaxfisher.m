% [R,E,J] = Bim_segmaxfisher(I,p)
%
% Toolbox: Balu
%  Segmentation of an object with homogeneous background. The idea is to
%  find a linear transformation from RGB to grayscale so that the Fisher
%  discriminant is maximum.
%
%  I: input image
%  p: threshold (default: p=-0.05) with p between -1 and 1.
%      A positive value is used to dilate the segmentation,
%      the negative to erode.
%  R: binary image of the object
%  E: binary image of the edge of the object
%  J: high contrast image of I.
%
%
%  Example:
%     I = imread('testimg1.jpg');
%     R = Bim_segmaxfisher(I);
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
% D.Mery, PUC-DCC, Apr. 2011
% http://dmery.ing.puc.cl

function [R,E,J] = Bim_segmaxfisher(I,p)
if size(I,3)==3
    J = Brgb2hci(double(I)/256);
else
    J = Bim_maxmin(I);
    n = fix(size(J,1)/4);
    if (mean2(J(1:n,1:n)) > 0.4)
        J = 1 - J;
    end
end
t = graythresh(J);
if (~exist('p','var'))
    p = -0.05;
end
[R,E] = Bim_morphoreg(J,t+p);
end


% J = Bim_rgb2hci(RGB)
%
%    RGB: color image
%    J  : high contrast image

function J = Brgb2hci(RGB)

RGB = double(RGB);
if (size(RGB,3)==1)
    I = RGB;
else
    RGB64 = imresize(RGB,[64 64]);
    k = fminsearch(@Bstd2,[1 1],[],RGB64);
    I = k(1)*RGB(:,:,1) + k(2)*RGB(:,:,2) + RGB(:,:,3);
end
J = Bim_maxmin(I);
n = fix(size(J,1)/4);
if (mean2(J(1:n,1:n)) > 0.4)
    J = 1 - J;
end
end


% s = Bstd2(k,RGB)
%
% Toolbox: Balu
%    Fisher discriminant of bimodal image I where
%    I = k(1)*R+k(2)*G+B (R = RGB(:,:,1), G=RGB(:,:,2), B = RGB(:,:,3)

function s = Bstd2(k,RGB)


I = k(1)*RGB(:,:,1) + k(2)*RGB(:,:,2) + RGB(:,:,3);
I = Bim_maxmin(I);
g = graythresh(I);
i1 = find(I>g);
i2 = find(I<=g);
m1 = mean(I(i1));
m2 = mean(I(i2));
v1 = var(I(i1));
v2 = var(I(i2));
s = -(m1-m2)/(v1+v2);

end