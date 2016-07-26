% [f1,d1,f2,d2,scores] = Bmv_matchSIFT(I1,I2,method,show)
%
% Toolbox: Balu
%
%    Matching points between two images I1 and I2 using SIFT points.
%
%    This function requires VLFeat Toolbox from (www.vlfeat.org).
%
%    method = 1 is for vl_ubcmatch method
%    method = 2 is for vl_ubcmatch plus RANSAC with homography
%               (modified from sift_mosaic.m by Andrea Vedaldi)
%
%    Example:
%       I1 = imread('testimg5.jpg');            % Image 1
%       I2 = imread('testimg6.jpg');            % Image 2
%       Bmv_matchSIFT(I1,I2,1,1)                % Method 1
%
%       Bmv_matchSIFT(I1,I2,2,1)                % Method 2
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [f1,d1,f2,d2,scores] = Bmv_matchSIFT(I1,I2,method,show)

if ~exist('show','var')
    show=0;
end

if show
    figure(1)
    clf
    imshow(I1)
    hold on
    figure(2)
    clf
    imshow(I2)
    hold on
    drawnow
end

if size(I1,3)==3
    im1 = single(rgb2gray(I1));
else
    im1 = single(I1);
end
if size(I2,3)==3
    im2 = single(rgb2gray(I2)) ;
else
    im2 = single(I2);
end

[f1, d1] = vl_sift(im1) ;
[f2, d2] = vl_sift(im2) ;
[matches, scores] = vl_ubcmatch(d1, d2) ;
[ii,jj] = sort(scores);
matches = matches(:,jj);
scores = scores(:,jj);

if method == 1


    if show
        figure(1)
        n = min([100 size(matches,2)]);
        disp(sprintf('%d best matching points',n'))
        for i=1:n
            plot(f1(1,matches(1,i)),f1(2,matches(1,i)),'*')
            text(f1(1,matches(1,i)),f1(2,matches(1,i)),num2str(i))
        end
        figure(2)
        for i=1:n
            plot(f2(1,matches(2,i)),f2(2,matches(2,i)),'*')
            text(f2(1,matches(2,i)),f2(2,matches(2,i)),num2str(i))
        end
    end


else
    %function [matches,scores] = Bsiftmatch2(im1,im2,show)
    % SIFT_MOSAIC Demonstrates matching two images using SIFT and RANSAC
    % Modified from sift_mosaic.m by Andrea Vedaldi
    % --------------------------------------------------------------------
    %                                                         SIFT matches
    % --------------------------------------------------------------------

    %     im1 = I1;
    %     im2 = I2;
    %
    %     if size(im1,3)==3
    %         im1 = rgb2gray(im1);
    %     end
    %     if size(im2,3)==3
    %         im2 = rgb2gray(im2);
    %     end
    %     [f1,d1] = vl_sift(im2single(im1)) ;
    %     [f2,d2] = vl_sift(im2single(im2)) ;
    %
    %     [matches,scores] = vl_ubcmatch(d1,d2) ;

    numMatches = size(matches,2) ;

    X1 = f1(1:2,matches(1,:)) ; X1(3,:) = 1 ;
    X2 = f2(1:2,matches(2,:)) ; X2(3,:) = 1 ;

    % --------------------------------------------------------------------
    %                                         RANSAC with homography model
    % --------------------------------------------------------------------

    clear H score ok ;
    for t = 1:100
        % estimate homograpyh
        subset = vl_colsubset(1:numMatches, 4) ;
        A = [] ;
        for i = subset
            A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
        end
        [U,S,V] = svd(A) ;
        H{t} = reshape(V(:,9),3,3) ;

        % score homography
        X2_ = H{t} * X1 ;
        du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
        dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
        ok{t} = (du.*du + dv.*dv) < 6*6 ;
        score(t) = sum(ok{t}) ;
    end

    [score, best] = max(score) ;
    H = H{best} ;
    ok = ok{best} ;


    % % --------------------------------------------------------------------
    % %                                                  Optional refinement
    % % --------------------------------------------------------------------
    %
    %     function err = residual(H)
    %         u = H(1) * X1(1,ok) + H(4) * X1(2,ok) + H(7) ;
    %         v = H(2) * X1(1,ok) + H(5) * X1(2,ok) + H(8) ;
    %         d = H(3) * X1(1,ok) + H(6) * X1(2,ok) + 1 ;
    %         du = X2(1,ok) - u ./ d ;
    %         dv = X2(2,ok) - v ./ d ;
    %         err = sum(du.*du + dv.*dv) ;
    %     end
    %
    % if exist('fminsearch') == 2
    %     H = H / H(3,3) ;
    %     opts = optimset('Display', 'none', 'TolFun', 1e-8, 'TolX', 1e-8) ;
    %     H(1:8) = fminsearch(@residual, H(1:8)', opts) ;
    % else
    %     warning('Refinement disabled as fminsearch was not found.') ;
    % end

    % --------------------------------------------------------------------
    %                                                         Show matches
    % --------------------------------------------------------------------

    if show
        figure(3) ; clf ;
        subplot(2,1,1) ;
        imagesc([im1 im2]) ;
        o = size(im1,2) ;
        line([f1(1,matches(1,:));f2(1,matches(2,:))+o], ...
            [f1(2,matches(1,:));f2(2,matches(2,:))]) ;
        title(sprintf('%d tentative matches', numMatches)) ;
        axis image off ;

        subplot(2,1,2) ;
        imagesc([im1 im2]) ;
        o = size(im1,2) ;
        line([f1(1,matches(1,ok));f2(1,matches(2,ok))+o], ...
            [f1(2,matches(1,ok));f2(2,matches(2,ok))]) ;
        title(sprintf('%d (%.2f%%) inliner matches out of %d', ...
            sum(ok), ...
            100*sum(ok)/numMatches, ...
            numMatches)) ;
        axis image off ;

        drawnow ;
    end

    matches = matches(:,ok);
    scores  = scores(:,ok);
end
f1 = f1(:,matches(1,:));
d1 = d1(:,matches(1,:));
f2 = f2(:,matches(2,:));
d2 = d2(:,matches(2,:));
