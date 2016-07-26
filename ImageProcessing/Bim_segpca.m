% [R,E,J] = Bim_segpca(I,p)
%
% Toolbox: Balu
%  Segmentation of an object with homogeneous background using the first
%  principal component of I.
%
%  I: input image
%  p: threshold (default: p=-0.05) with p between -1 and 1.
%      A positive value is used to dilate the segmentation,
%      the negative to erode.
%  R: binary image of the object
%  E: binary image of the edge of the object
%  J: high contrast image of I.
%
%  Example:
%     I = imread('testimg1.jpg');
%     R = Bim_segpca(I);
%     figure(1)
%     imshow(I); title('test image')
%     figure(2)
%     imshow(R); title('segmented image')
%
%     Repeat this examples for images testimg2, testimg3 and testimg4. Last
%     test image requires R = Bim_segbalu(I,-0.1) for better results.
%
% See also Bim_segmowgli, Bim_segotsu, Bim_segkmeans, Bio_segshow,
% Bim_segbalu.
%
% D.Mery, PUC-DCC, May 2011
% http://dmery.ing.puc.cl

function [R,E,J] = Bim_segpca(I,p)


I   = double(I);



% PCA of whole image (time computing is very high)
% I1   = I(:,:,1);
% I2   = I(:,:,2);
% I3   = I(:,:,3);
% X    = [Bft_norm(I1(:),0) Bft_norm(I2(:),0) Bft_norm(I3(:),0)];
% Y    = Bft_pca(double(X),1);
% J    = zeros(size(I1));
% J(:) = Bft_norm(Y,0);

if size(I,3)==3

    Ir   = imresize(I,[64 64]);
    I1   = Ir(:,:,1);
    I2   = Ir(:,:,2);
    I3   = Ir(:,:,3);
    X    = [Bft_norm(I1(:),0) Bft_norm(I2(:),0) Bft_norm(I3(:),0)];

    [Y,lambda,A] = Bft_pca(X,1);

    I1t  = I(:,:,1);
    I2t  = I(:,:,2);
    I3t  = I(:,:,3);

    I1n  = I1t(:)-mean(I1(:));
    I2n  = I2t(:)-mean(I2(:));
    I3n  = I3t(:)-mean(I3(:));
    X0   = [I1n I2n I3n];
    B    = A(:,1);
    Yn   = X0*B; % Y = X0*A; Y = Y(:,1:m); % the first m components

    J = zeros(size(I1t));
    J(:) = Bft_norm(Yn,0);
else
    J = Bim_maxmin(I);

end
n = fix(size(J,1)/4);
if (mean2(J(1:n,1:n)) > 0.4)
    J = 1 - J;
end

t = graythresh(J);
if (~exist('p','var'))
    p = -0.05;
end
[R,E] = Bim_morphoreg(J,t+p);