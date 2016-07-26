% Y = Bim_color2bwreg(I,R)
%
% Toolbox: Balu
%   Convert only a region of a color image into a grayscale image.
%
%   I is a color image, R is a binary image, I and R have the same size.
%   This function converts to grayscale those pixels of I where R is equal 
%   to '1'. The rest of the pixels keep the same original color. The output 
%   is Y and it has the same size than I.
%
%  Example:
%     I = imread('testimg9.jpg');
%     R = not(imdilate(and(I(:,:,1)<70,I(:,:,2)<70),ones(5,5)));
%     J = Bim_color2bwreg(I,R);
%     figure(1)
%     imshow(I); title('original image')
%     figure(2)
%     imshow(R); title('mask')
%     figure(3)
%     imshow(J); title('converted image')
%     disp('The palm only is b&w, the sky and the clouds are colored.')
%
% (c) D.Mery, PUC-DCC, 2014
% http://dmery.ing.puc.cl




function Y = Bim_color2bwreg(I,R)

I        = uint8(I);
IR       = I(:,:,1);
IG       = I(:,:,2);
IB       = I(:,:,3);
K        = rgb2gray(I);
z        = R==0;
IR(z)    = K(z);
IG(z)    = K(z);
IB(z)    = K(z);
Y        = uint8(zeros(size(I)));
Y(:,:,1) = IR;
Y(:,:,2) = IG;
Y(:,:,3) = IB;
