% [Ibs,H] = Bmv_homographySIFT(Ia,Ib,Ra,Rb,show)
%
% Toolbox: Balu
%
%    Homography between images Ia and Ib using RANSAC of SIFT points
%    Ibs is transformed image Ib.
%    Ra and Rb are binary images that indicate where the keypoints are
%    valid, ie, SIFT descriptors in pixels of Ia (or Ib) where Ra = 0 (or
%    Rb = 0) will not be considered.
%    size(Ia) = size(Ibs)
%    If ma and mb are the homogeneus coordinates of points in images Ia and Ib:
%       ma = [xa ya 1]'; mb = [xb yb 1]';
%       then  mb = H*ma.
%    show = 1 displays results.
%
%   Example:
%      Ia = rgb2gray(imread('testimg4.jpg'));           % Image Ia
%      Ra = Bim_segbalu(Ia);                            % segmentation of Ia
%      H  = [1.9 0.01 0.01;0.002 1.95 0.001;0.001 0 1]  % Original H
%      Ib = Bmv_projective2D(Ia,H,[1000 750],1);        % Image Ib
%      Rb = Bim_segbalu(Ib);                            % segmentation of Ib
%      [Ibs,Hs] = Bmv_homographyRSIFT(Ia,Ib,Ra,Rb,1);   % Homography
%      Hs = Hs/Hs(3,3)                                  % estimated matrix H
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [Ibs,H] = Bmv_homographyRSIFT(Ia,Ib,Ra,Rb,show)

if ~exist('show','var')
    show = 0;
end


if isempty(Ra)
    Ra = ones(size(Ia));
end
if isempty(Rb)
    Rb = ones(size(Ib));
end


if show
    figure(4)
    clf
    imshow(Ia,[])
    hold on
    figure(5)
    clf
    imshow(Ib,[])
    hold on
    drawnow
    disp('searching SIFT points...');
end
% nsize      = options.nsize ;
% binsize    = options.binsize ;
% magnif     = options.magnif ;
% sfactor    = sqrt((binsize/magnif)^2 - .25);
% 
% Iaf     = vl_imsmooth(single(Ia), sfactor);
% [fa0, da0] = vl_dsift(Iaf, 'size', binsize,'step',nsize,'fast');
% [fa,da] = andsift(fa0,da0,Ra);
% 
% Ibf     = vl_imsmooth(single(Ib), sfactor);
% [fb0, db0] = vl_dsift(Ibf, 'size', binsize,'step',nsize,'fast');
% [fb,db] = andsift(fb0,db0,Rb);

Ia = single(Ia) ;
Ib = single(Ib) ;
[fa0, da0] = vl_sift(Ia) ;
% [fa0, da0] = vl_dsift(Ia) ;
[fa,da] = andsift(fa0,da0,Ra);
[fb0, db0] = vl_sift(Ib) ;
% [fb0, db0] = vl_dsift(Ib) ;
[fb,db] = andsift(fb0,db0,Rb);
[matches, scores] = vl_ubcmatch(da, db) ;

[ii,jj] = sort(scores); %#ok<ASGLU>
matches = matches(:,jj);

n = min([100 size(matches,2)]);

xa = fa(2,matches(1,1:n))';
ya = fa(1,matches(1,1:n))';

xb = fb(2,matches(2,1:n))';
yb = fb(1,matches(2,1:n))';

ma = [xa'; ya'; ones(1,n)];
mb = [xb'; yb'; ones(1,n)];

if show
    disp('computing RANSAC homography...');
end

H = Bmv_homographyRANSAC(ma,mb);


if show
    
    figure(4)
    fprintf('%d most similar points...\n',n')
    for i=1:n
        plot(fa(1,matches(1,i)),fa(2,matches(1,i)),'*')
        text(fa(1,matches(1,i)),fa(2,matches(1,i)),num2str(i))
    end
    figure(5)
    for i=1:n
        plot(fb(1,matches(2,i)),fb(2,matches(2,i)),'*')
        text(fb(1,matches(2,i)),fb(2,matches(2,i)),num2str(i))
    end
end

Ibs = Bmv_projective2D(Ib,H,size(Ia),show);

ni = sum(Ibs(:)==0)/length(Ibs(:));

if show
    figure(3)
    
    subplot(2,2,1)
    imshow(Ia,[])
    xlabel('Image Ia')
    
    subplot(2,2,2)
    imshow(Ib,[])
    xlabel('Image Ib')
    
    subplot(2,2,3)
    imshow(Ibs,[])
    xlabel('Ibs: Transformed Image Ib')
    
    subplot(2,2,4)
    imshow(double(Ia)+Ibs,[])
    xlabel('(Ia + Ibs)/2')
end

if ni>0.5
    error('Bmv_homographyRSIFT: bad estimation of matrix H.')
end

