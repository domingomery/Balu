% function Y = LUT(X,T,show)
%
% Toolbox: Balu
%    Look Up Table transformation for grayvalue images.
%
%    Input data:
%       X grayvalue image.
%       T look uo table
%       show display results
%
%    Output:
%       Y transformed image
%
% Ejemplo:
%    load clown
%    T = 256*ones(256,1);
%    T(1:81) = (1:81).*(1:81)/81^2*255;; % look up table
%    Y = Bim_LUT(X,T,1); 
%
% (c) D.Mery, PUC-DCC, 2009
% http://dmery.ing.puc.cl

function Y = Bim_LUT(X,T,show)

if ~exist('show','var')
    show = 0;
end



X = fix(double(X)+0.5);

ii = find(X>256);
if not(isempty(ii))
    X(ii)=256*ones(length(ii),1);
end

ii = find(X<1);
if not(isempty(ii))
    X(ii)=ones(length(ii),1);
end

Y = uint8(T(X));

if not(exist('show','var'))
    show = 0;
end

if (show)
    close all
    figure
    imshow(X,gray(256));title('original image');
    figure
    imshow(Y,gray(256));title('transformaded image');

    figure
    plot(T)
    axis([0 256 0 256])
    title('LUT function');
    xlabel('input');
    ylabel('output');

    figure
    hist(double(X(:)),0:255);title('original histogram');
    ax = axis;
    ax(1:2) = [0 256];
    axis(ax)

    figure
    hist(double(Y(:)),0:255);title('transformed histogram');
    ax = axis;
    ax(1:2) = [0 256];
    axis(ax)
    
    figure
    subplot(3,3,1)
    plot(T)
    axis([0 256 0 256])
    title('LUT function');
    xlabel('input');
    ylabel('output');

    subplot(3,3,4)
    y = hist(double(X(:)),0:255);
    bar(y);
    title('original histogram');
    ax = axis;
    ax(1:2) = [0 256];
    axis(ax)


    subplot(3,3,2)
    y = hist(double(Y(:)),0:255);title('transformed histogram');
    barh(y)
    ax = axis;
    ax(3:4) = [0 256];
    axis(ax)
    
    subplot(3,3,3)
    imshow(Y,gray(256))
    title('transformed image');

    subplot(3,3,7)
    imshow(X,gray(256))
    title('original image');
end