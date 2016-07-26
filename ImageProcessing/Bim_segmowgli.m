%  [F,m] = Bsegmowgli(J,R,Amin,sig)
% 
% Toolbox: Balu
%  Segmentation of regions in image J using LoG edge detection.
%  R   : binary image of same size of J that indicates the piexels where
%        the segmentation will be performed. Default R = ones(size(J));
%  Amin: minimum area of the segmented details.
%  sig : sigma of LoG edge detector.
%  F   : labeled image of the segmentation.
%  m   : numbers of segmented regions.
% 
%  Example 1:
%     I = imread('rice.png');
%     figure(1);imshow(I);title('test image');
%     [F,m] = Bim_segmowgli(I,[],40,1.5);
%     figure(2);imshow(F,[]);title('segmented image');     
% 
%
%  Example 2:
%     I = imread('testimg4.jpg');
%     figure(1);imshow(I);title('test image');
%     R = Bim_segbalu(I,-0.1);
%     figure(2);imshow(R);title('segmented object');
%     G = I(:,:,2);
%     [F,m] = Bim_segmowgli(G,R,30,2);
%     figure(3);imshow(F,[]);title('segmented regions')
%
%     Another way to display the results:
%        Bio_edgeview(I,F>0)
%
% See also Bsegbalu.
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl
%

function [F,m] = Bim_segmowgli(J,R,Amin,sig)
if (not(exist('R','var')))
    R = ones(size(J));
end
if isempty(R)
   R = ones(size(J));
end 
if (not(exist('sig','var')))
    sig = 2;
end
if (not(exist('Amin','var')))
    Amin = 20;
end
[N,M,P] = size(J);
if (P==3)
    J = rgb2gray(J);
end
se = strel('disk',3);
Re = imdilate(R,se);
E  = bwperim(R,4);
L  = (edge(J,'log',1e-10,sig) & Re) | E;
W  = bwareaopen(L,Amin);
F  = bwlabel(not(W),4);
m  = max(max(F));
for i=1:m
    ii = find(F==i);
    if (length(ii)<Amin)||(length(ii)>N*M/6)
        W(ii) = 1;
    end
end
F  = bwlabel(not(W),4);
m  = max(F(:));
fprintf('%d segmented regions.\n',m);
