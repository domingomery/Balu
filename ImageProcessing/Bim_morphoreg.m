% [R,E] = Bmorphoreg(J,t);
% [R,E] = Bmorphoreg(Ro);
%
% Toolbox: Balu
%    Morphology operations of binary image J>t (or Ro): remove isolate 
%    pixels and fill holes.
%    R: binary image of the region
%    E: binary image of the edge
%
%  Example:
%     I = imread('testimg2.jpg');
%     figure(1);imshow(I)
%     J = rgb2gray(I);
%     Ro = Bim_segotsu(J);
%     figure(2);imshow(Ro)
%     [R,E] = Bim_morphoreg(Ro);
%     figure(3);imshow(R)
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl
%

function [R,E] = Bim_morphoreg(J,t)

if ~exist('t','var')
    Ro = J;
else
    Ro = J>t;
end



A = bwareaopen(Ro,fix(length(Ro(:))/100));
C = imclose(double(A),strel('disk',7));
R = bwfill(C,'holes',8);
E = bwperim(R,4);