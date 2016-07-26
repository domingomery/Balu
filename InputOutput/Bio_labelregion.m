% [d,D] = Bio_labelregion(I,L,c)
%
% Toolbox: Balu
%    User interface to label regions of an image.
%
%    I is the original image (color or grayvalue).
%    L is a labeled image that indicates the segmented regions of I.
%    c is the maximal number of classes.
%    d(i) will be the class number of region i.
%    D is a binary image with the corresponding labels.
%
% Example:
%   I = imread('rice.png');                          
%   I = I(1:70,1:70);                              % input image
%   [L,m] = Bim_segmowgli(I,ones(size(I)),40,1.5); % segmented image
%   [d,D] = Bio_labelregion(I,L,3);
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl
%


function [d,D] = Bio_labelregion(I,L,c)


if size(I,3)==3
    J = rgb2gray(I);
else
    J = I;
end



close all


cmap = [0 0 1
    0 1 0
    1 0 0
    0 1 1
    1 0 1
    1 1 0
    0 0 0
    1 1 1];

cmap = [cmap;cmap;cmap];

figure(2)
Icol = [];
for i=1:c
    Icol = [Icol;ones(20,20)*i];
end
imshow(Icol,cmap)
hold on
for i=1:c
    text(5,i*20-10,num2str(i))
end
title('Label color')
colorstr = 'bgrcmykwbgrcmykwbgrcmykw';

warning off
figure(1)
subplot(1,2,2)
imshow(I,[])
title('Original image')
hold on

n = max(max(L));
d = zeros(n,1);

i = 1;
D = zeros(size(J));


while(i<=n)

    R = zeros(size(J));
    kk = find(L==i);
    R(kk) = 1;
    E = bwperim(R,4);

    %subplot(2,2,3)
    %Bio_edgeview(J,E);
    %title('Class label?')

    %figure(3)
    subplot(1,2,1)
    Bio_edgeview(J,R,[1 1 0]);
    title('Class label of yellow region?')


    ok = 0;

    while(not(ok))
        r = input(sprintf('Region %3d/%3d: Class label? (1... %d) or -1 to correct previous: ',i,n,c));

        % r = input(['Class label? (1...' num2str(c) ') or -1 to correct previous: ']);
        r = round(r);
        if not(isempty(r))
            if (r==-1)
                i = max(i-2,0);
                ok = 1;
                disp('Correction: enter new choise...')
            end
            if (r>=1) && (r<=c)
                d(i) = r;
                D(kk) = r;
                ok = 1;
                [ii,jj] = find(R==1);
                x1 = max(jj);
                x2 = min(jj);
                y1 = max(ii);
                y2 = min(ii);
                subplot(1,2,2)
                plot([x1 x1 x2 x2 x1],[y1 y2 y2 y1 y1],colorstr(r));
                title('Labeled regions')
            end
        end
        if (ok==0)
            beep
        end
    end
    i = i + 1;
end

Dc = D+1;
map = [0 0 0;cmap];
subplot(1,2,2)
imshow(Dc,map);
subplot(1,2,1)
Bio_edgeview(J,zeros(size(J)),[1 1 0]);
title('Original image')