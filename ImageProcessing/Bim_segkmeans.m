% [R,J] = Bim_segkmeans(I,k,r,show)
%
% Toolbox: Balu
%
%    Segmentation of color images using kmeans.
%
% I: input image
% k: number of clusters
% r: resize image
% show: 1 means intermediate results will be displayed
%
% Example 1:
%    I = imread('testimg1.jpg');
%    [R,J] = Bim_segkmeans(I,2,1,1);
% 
% Example 2:
%   I = imread('testimg9.jpg');
%   [R,J] = Bim_segkmeans(I,3,1,1);
%
%  See also Bim_segbalu, Bim_segotsu.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [R,J] = Bim_segkmeans(I,k,r,show)

I = double(I);
if (exist('r','var'))
    I = imresize(I,r);
end

if show
    close all
    imshow(I/256)
    title('Original Image')
    drawnow
end
[N,M,P] = size(I);

if P~=3
    error('Bim_segkmeans works only with 3 channel images like RGB.');
end


R = zeros(N,M);


xr = I(:,:,1);
xg = I(:,:,2);
xb = I(:,:,3);

X = [xr(:) xg(:) xb(:)];


ds = Bct_kmeans(X,k);
R(:) = ds;
R = uint8(R);

xr = I(:,:,1);
xg = I(:,:,2);
xb = I(:,:,3);


yr = zeros(size(xr));
yg = yr;
yb = yr;

for i=1:k
    ii = find(ds==i);

    if not(isempty(ii))
        yr(ii) = mean(xr(ii));
        yg(ii) = mean(xg(ii));
        yb(ii) = mean(xb(ii));
        if show
            K = zeros(size(I));
            kr = K(:,:,1);
            kg = K(:,:,2);
            kb = K(:,:,3);
            kr(ii) = xr(ii);
            kg(ii) = xg(ii);
            kb(ii) = xb(ii);
            K(:,:,1) = kr;
            K(:,:,2) = kg;
            K(:,:,3) = kb;
            figure(2)
            imshow(uint8(K))
            title(sprintf('Cluster %d',k))
            enterpause
        end
    end
end
J = zeros(size(I));
J(:,:,1) = yr;
J(:,:,2) = yg;
J(:,:,3) = yb;
J = uint8(fix(J));

if show
    figure(2)
    imshow(J)
    title('Segmented Image');
end
