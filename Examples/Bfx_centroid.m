% [X,Xn] = Bfx_centroid(R,options)
%
% Toolbox: Balu
%
%    Centroid of a region.
%
%    options.show    = 1 display mesagges.
%
%      X(1) is centroid-i, X(2) is centroid-j
%      Xn is the list of the n feature names.
%
%    Example (Centroid of a region)
%      I = imread('testimg1.jpg');     % input image
%      R = Bim_segbalu(I);             % segmentation
%      imshow(R);
%      options.show = 1;
%      [X,Xn] = Bfx_centroid(R,options);
%      Bio_printfeatures(X,Xn)
%
%   See also Bfx_basicgeo, Bfx_gupta, Bfx_fitellipse, Bfx_flusser, 
%   Bfx_hugeo.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl

function [X,Xn] = Bfx_centroid(R,options)

if options.show == 1
    disp('--- extracting centroid...');
end
[Ireg,Jreg] = find(R==1);           % pixels in the region
ic = mean(Ireg);
jc = mean(Jreg);
X  = [ic jc];
Xn = [ 'Centroid i              '
       'Centroid j              ']; % 24 characters per name
if options.show == 1
    clf
    imshow(R)
    hold on
    plot(X(2),X(1),'rx')
    enterpause
end
