% Bio_edgeview(B,E,c,g)
%
% Toolbox: Balu
%   Display gray or color image I overimposed by color pixels determined
%   by binary image E. Useful to display the edges of an image.
%   Variable c is the color vector [r g b] indicating the color to be displayed 
%  (default: c = [1 0 0], i.e., red)
%   Variable g is the number of pixels of the edge lines, default g = 1
%
% Example to display a red edge of a food: 
%    I = imread('testimg1.jpg');           % Input image
%    [R,E] = Bim_segbalu(I);               % Segmentation
%    Bio_edgeview(I,bwperim(E),[0 1 0],3)  % perimeter with 3 pixels will be
%                                          % displayed ingreen ([0 1 0] for [R G B])
%
% D.Mery, PUC-DCC, Apr. 2008-2019
% http://dmery.ing.puc.cl
%

function Bio_edgeview(B,E,cc,g)

if not(exist('cc','var'))
    cc = [1 0 0];
end

if not(exist('g','var'))
    g = 1;
end


B = double(B);
if max(B(:))>1
    B = B/256;
end

if (size(B,3)==1)
    [N,M] = size(B);
    J = zeros(N,M,3);
    J(:,:,1) = B;
    J(:,:,2) = B;
    J(:,:,3) = B;
    B = J;
end

B1 = B(:,:,1);
B2 = B(:,:,2);
B3 = B(:,:,3);

Z = B1==0;
Z = and(Z,B2==0);
Z = and(Z,B3==0);
ii = find(Z==1);
if not(isempty(ii))
    B1(ii) = 1/256;
    B2(ii) = 1/256;
    B3(ii) = 1/256;
end
warning off
E = imdilate(E,ones(g,g));
ii       = find(E==1);
B1(ii)   = cc(1)*255;
B2(ii)   = cc(2)*255;
B3(ii)   = cc(3)*255;
Y        = double(B);
Y(:,:,1) = B1;
Y(:,:,2) = B2;
Y(:,:,3) = B3;
imshow(uint8(Y*256))
drawnow
warning on
