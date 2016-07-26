% J = Bimgd2(I)
%
% Toolbox: Balu
%    Second derivative of image X.
%
%    Input data:
%       I grayvalue image.
%
%    Output:
%       J = conv2(I,[0 1 0;1 -4 1;0 1 0],'same'); 
%
%    Example:
%       X = imread('testimg2.jpg');
%       I = rgb2gray(X);
%       figure(1)
%       imshow(I); title('original image')
%       J = Bimgd2(I);
%       figure(2)
%       imshow(J,[]); title('laplacian')
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function Y = Bim_d2(X)
Y = conv2(X,[0 1 0;1 -4 1;0 1 0],'same');
