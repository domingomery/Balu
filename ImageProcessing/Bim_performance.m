% [precision,recall] = Bim_performance(Xideal,Xtest,thfp,thtp,show)
%
% Toolbox: Balu
%
%    Precision and recall rates using binary images Xideal (for the ideal
%    detection) and Xtest (for the real detection). 
%
%    The computation of the False Positives are as follows: all pixels from 
%    Xtest that have a minimal distance from all pixels of Xideal larger 
%    than thfp are marked. These pixels are dilated using a disk mask with
%    radius thfp/2. The number of these pixels are the number of FP.
%
%    The computation of True Positives are as follows:  all pixels from 
%    Xideal that have a minimal distance from all pixels of Xtest less 
%    than thtp are marked. The number of these pixels are the number of TP.
%
%    The performance rates are defined as:
%      Precision = TP/(TP+FP)
%      Recall    = TP/P     
%
% Example:
%    X      = imread('testimg7.bmp'); % grayscale image
%    figure(1)
%    imshow(X,[]);
%    Xideal = imread('testimg8.bmp'); % ideal detection
%    X1     = Bim_segotsu(X);
%    X2     = imerode(X1,ones(35,35));
%    Xtest  = and(X2,X<120);
%    [precision,recall] = Bim_performance(Xideal,Xtest,10,7,1);%
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [precision,recall] = Bim_performance(Xideal,Xtest,thfp,thtp,show)

[xi,yi] = find(Xideal==1);
[xj,yj] = find(Xtest==1);
ci      = [xi'; yi']; %ideal detection
cj      = [xj'; yj']; % real detection
ki      = vl_kdtreebuild(ci);
[ji,Di] = vl_kdtreequery(ki,ci,cj,'NumNeighbors',1);
di      = sqrt(Di);

Npos    = length(xi);
ifp     = find(di>thfp);

ZFP     = zeros(size(Xtest));
k       = sub2ind(size(ZFP),xj(ifp),yj(ifp));
ZFP(k)  = 1;

se      = strel('disk',round(thfp/2));
ZFPd    = imdilate(ZFP,se);

FP      = sum(ZFPd(:));

kj      = vl_kdtreebuild(cj);
[jj,Dj] = vl_kdtreequery(kj,cj,ci,'NumNeighbors',1);
dj      = sqrt(Dj);

itp     = find(dj<thtp);

TP      = length(itp);
ZTP     = uint8(Xideal);
k       = sub2ind(size(ZTP),xi(itp),yi(itp));
ZTP(k)  = 2; %*ones(length(k),1);

recall    = TP/Npos;
precision = TP/(TP+FP);
if show
    figure(2)
    clf
    imshow(Xideal)
    hold on
    plot(yj,xj,'r.');
    title('Ideal (white) and real detection (red)')
    figure(3)
    imshow(ZFPd,[]);
    title('False detections (FP:white)')
    figure(4)
    imshow(ZTP,[]);
    title('True positives (TP:white) and Missdetections (FN:gray)')
    fprintf('Precision = TP/(TP+FP) = %7.4f\n',precision)
    fprintf('Recall    = TP/P       = %7.4f\n',recall)
end






