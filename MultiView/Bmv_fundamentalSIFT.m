% F = Bmv_fundamentalSIFT(I1,I2)
%
% Toolbox: Balu
%
%    Estimation of Fundamental Matrix from two images using SIFT points.
%
%    This function requires VLFeat Toolbox from (www.vlfeat.org).
%
%    I1 and I2 are the stereo images.
%    F is the Fundamental matrix.
%
%    Example:
%       I1 = imread('testimg5.jpg');            % Image 1
%       figure(1);imshow(I1); hold on
%       I2 = imread('testimg6.jpg');            % Image 2
%       figure(2);imshow(I2); hold on
%       F  = Bmv_fundamentalSIFT(I1,I2);        % F matrix estimation
%       while(1)
%          figure(1);
%          disp('click a point in Figure 1...') % click
%          p = vl_click; m1 = [p(2) p(1) 1]';
%          plot(p(1),p(2),'g+')
%          figure(2)
%          Bmv_epiplot(F,m1)
%       end
%
%    See also Bmv_epiplot, Bmv_fundamentalRANSAC, Bmv_fundamentalSVD. 
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function F = Bmv_fundamentalSIFT(I1,I2)

if ~exist('vl_sift','file')
    error('Bmv_fundamentalSIFT: This function requires the VLFeat Toolbox.');
end

if size(I1,3)==3;
    I1 = rgb2gray(I1);
end

if size(I2,3)==3;
    I2 = rgb2gray(I2);
end

im1 = single(I1);
im2 = single(I2);

% extracting SIFT points...
[f1, d1] = vl_sift(im1) ;
[f2, d2] = vl_sift(im2) ;

frap = f1';
fraq = f2';

% searching matching points...
[matches, scores] = vl_ubcmatch(d1, d2);

ii       = find(scores<90000);
xp       = frap(matches(1,ii),[2 1]);
xq       = fraq(matches(2,ii),[2 1]);

% computing fundamental matrix using RANSAC...
n = size(xp,1);
m1 = [xp';ones(1,n)];
m2 = [xq';ones(1,n)];
F = Bmv_fundamentalRANSAC(m1,m2);

