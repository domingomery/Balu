% K = Bim_rgb2pca(I)
%
% Toolbox: Balu
%    Converts a RGB image into a new color image where the channel i is
%    the principal componet i of image I, for i=1,2,3.
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl


function K = Bim_rgb2pca(I)

I    = double(I);
X1   = I(:,:,1);
X2   = I(:,:,2);
X3   = I(:,:,3);
X    = [X1(:) X2(:) X3(:)];
Y    = Bft_pca(X,3);
Y1    = zeros(size(X1));
Y2    = zeros(size(X1));
Y3    = zeros(size(X1));
Y1(:) = Y(:,1);
Y2(:) = Y(:,2);
Y3(:) = Y(:,3);
K = zeros(size(I));
K(:,:,1) = Y3;
K(:,:,2) = Y2;
K(:,:,3) = Y1;
