function [xy,bs,pose,Icrop] = Bfr_posedetection(I,show)


% X. Zhu, D. Ramanan. "Face detection, pose estimation and landmark
% localization in the wild" Computer Vision and Pattern Recognition (CVPR)
% Providence, Rhode Island, June 2012.
%
% see domingomery/Matlab/face-release1.0-basic
%
% http://www.ics.uci.edu/~xzhu/face/
%
% This code was written by Domingo Mery from demo.m

xy = 0;
bs = 0;
pose = 1000;
Icrop = [];
if ~exist('show','var')
    show = 0;
end

% load and visualize model
% Pre-trained model with 146 parts. Works best for faces larger than 80*80
load face_p146_small.mat

% % Pre-trained model with 99 parts. Works best for faces larger than 150*150
% load face_p99.mat

% % Pre-trained model with 1050 parts. Give best performance on localization, but very slow
% load multipie_independent.mat


% 5 levels for each octave
model.interval = 5;
% set up the threshold
model.thresh = min(-0.65, model.thresh);

% define the mapping from view-specific mixture id to viewpoint
if length(model.components)==13
    posemap = 90:-15:-90;
elseif length(model.components)==18
    posemap = [90:-15:15 0 0 0 0 0 0 -15:-15:-90];
else
    error('Can not recognize this model');
end

if show>0
    clf; imagesc(I); axis image; axis off; drawnow;hold on
end

bs = detect(I, model, model.thresh);
bs = clipboxes(I, bs);
bs = nms_face(bs,0.3);

if ~isempty(bs)
    
    pose = posemap(bs(1).c);
    
    x = bs(1).xy(:,[1 3],:);y = bs(1).xy(:,[2 4]);x1 = min(x(:));x2 = max(x(:));y1 = min(y(:));y2 = max(y(:));
    
    dx = x2-x1;
    dy = y2-y1;
    
    x1 = round(x1+0.1*dx);
    x2 = round(x2-0.1*dx);
    
    y1 = round(y1-0.1*dy);
    y2 = round(y2-0.1*dy);
    
    
    
    xy = [x1 y1 x2 y2];
    Icrop = I(y1:y2,x1:x2,:);
    if show>0
        if show==2
            % show highest scoring one
            figure(1)
            showboxes(I, bs(1),posemap);
            title('Highest scoring detection');
            % figure,showboxes(I, bs,posemap),title('All detections above the threshold');
            % show all
        end
        plot([x1 x1 x2 x2 x1],[y1 y2 y2 y1 y1],'r');
        figure(2)
        imshow(Icrop);
    end
    
    
    
else
    bs = -1;
    pose = 1000;
end