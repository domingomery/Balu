% R = Bim_segotsu(I)
%
% Toolbox: Balu
%
%    Otsu segmentation of a grayvalue. This function requires Image
%    Processing Toolbox.
%
%    Input data:
%       I grayvalue image.
%
%    Output:
%       R: binary image.
%
%    Example: Training & Test together:
%       X = imread('testimg1.jpg');
%       figure(1)
%       imshow(X); title('original image')
%       R = Bim_segotsu(I);
%       figure(2)
%       imshow(R); title('segmented image')
%
%  See also Bim_segbalu, Bim_segkmeans.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [R,E,J] = Bim_segotsu(I,p)
if ~exist('p','var')
    p = 0;
end

Id = double(I);
if size(I,3)==3
    Id = rgb2gray(Id/256);
end
J = Bim_maxmin(Id);
n = fix(size(J,1)/4);
if (mean2(J(1:n,1:n)) > 0.4)
    J = 1 - J;
end

t = graythresh(J);
[R,E] = Bim_morphoreg(J,t+p);