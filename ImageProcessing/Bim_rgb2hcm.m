% J = Bim_rgb2hcm(RGB)
%
% Toolbox: Balu
%    Conversion RGB to high contrast image.
%
%    RGB: color image
%    J  : hcm image
%
%  See details in:
%  Mery, D.; Pedreschi, F. (2005): Segmentation of Colour Food Images using
%  a Robust Algorithm. Journal of Food Engineering 66(3): 353-360.
%
%  Example:
%     I = imread('testimg2.jpg');
%     J = Bim_rgb2hcm(I);
%     figure(1)
%     imshow(I); title('control image')
%     figure(2)
%     imshow(J); title('high contrast image')
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl


function J = Bim_rgb2hcm(RGB)

RGB = double(RGB);
if (size(RGB,3)==1)
    I = RGB;
else
    RGB64 = imresize(RGB,[64 64]);
    k = fminsearch(@Bstdmono,[1 1],[],RGB64);
    I = k(1)*RGB(:,:,1) + k(2)*RGB(:,:,2) + RGB(:,:,3);
end
J = I - min(I(:));
J = J/max(J(:));
[N,M] = size(J);
n = min([fix(size(J,1)/4) N]);
m = min([fix(size(J,2)/4) M]);

if (mean2(J(1:n,1:m)) > 0.4)
    J = 1 - J;
end
end


% s = Bstdmono(k,RGB)
%
% Toolbox: Balu
%    Standard deviation of normalized image I, where
%    I = k(1)*R+k(2)*G+B (R = RGB(:,:,1), G=RGB(:,:,2), B = RGB(:,:,3)
%    s: Standard deviation
%
%    This function is called by Brgb2hcm, Bsegbalu.
%
% See also Brgb2hcm, Bsegbalu.
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl
%

function s = Bstdmono(k,RGB)
I = k(1)*RGB(:,:,1) + k(2)*RGB(:,:,2) + RGB(:,:,3);
s = -std2(I)/(max(I(:))-min(I(:)));
end