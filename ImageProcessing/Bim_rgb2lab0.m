% Lab = Bim_rgb2lab0(RGB)
%
% Toolbox: Balu
%    Conversion RGB to L*a*b using formulas.
%
%    RGB: color image
%    Lab  : L*a*b*
%
%  Example:
%     I = imread('testimg2.jpg');
%     J = Bim_rgb2lab0(I);
%     figure(1)
%     imshow(I); title('control image')
%     figure(2); imshow(J(:,:,1),[]); title('L*')
%     figure(3); imshow(J(:,:,2),[]); title('a*')
%     figure(4); imshow(J(:,:,3),[]); title('b*')
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function Lab = Bim_rgb2lab0(RGB)
if max(RGB(:))>1
    RGB = double(RGB)/256;
end
XYZ = vl_rgb2xyz(RGB);
Lab = vl_xyz2lab(XYZ);