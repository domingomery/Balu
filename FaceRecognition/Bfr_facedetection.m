function [xy,Icrop] = Bfr_facedetection(I,show)

faceDetector = vision.CascadeObjectDetector('FrontalFaceCART');

bbox = step(faceDetector, I);

[N,M] = size(I);

s = size(bbox,1);

if s==0 % no detection
    xy = [];
    Icrop = [];
else
    if s>1 % it selects the largest detection
        [~,jj] = max(bbox(:,3));
        bbox = bbox(jj(1),:);
    end
    
    
    dx = bbox(3);
    dy = round(bbox(3)*11/9); % aspect 110 x 90
    
    d  = round((dy-bbox(4))/2);
    
    
    xy = bbox([1 1 2 2])+[0 dx -d dy-d];
    
    x1=xy(1);x2=xy(2);y1=max([1 xy(3)]);y2=min([xy(4) N]);
    xy(3) = y1; xy(4) = y2;
    Icrop = I(y1:y2,x1:x2,:);
    
    bbox = [bbox(1) max([1 bbox(2)-d]) dx dy];
    
    
    % Draw boxes around detected faces and display results
    if show==1
        figure(1),imshow(I);
        IFaces = insertObjectAnnotation(I, 'rectangle', bbox, 'Face');
        figure(2), imshow(IFaces), title('Detected Faces');
    end
end