% Y = Bim_equalization(X)
% Y = Bim_equalization(X,show)
%
% Toolbox: Balu
%    Enhancement of a grayvalue image forcing an uniform histogram.
%
%    Input data:
%       X grayvalue image.
%
%    Output:
%       Y: enhanced image so that histogram of Y is perfectlty uniformed
%       distributed. Y is uint8.
%
%    Example:
%       X = imread('tire.tif');
%       Y = Bim_equalization(X,1);
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl


function Y = Bim_equalization(X,show)

if ~exist('show','var')
    show = 0;
end


X     = uint8(round(Bim_lin(X)));

[N,M] = size(X);
Y     = uint8(zeros(N,M));
[i,j] = sort(X(:));

z     = uint8(zeros(N*M,1));
d     = fix(N*M/256+0.5);
for i=1:255;
    z((i-1)*d+1:i*d) = (i-1)*ones(d,1);
end
z(255*d+1:N*M) = 255*ones(N*M-255*d,1);

Y(j) = z;



if show
    figure(1)
    subplot(2,1,1);imshow(X); title('original image')
    subplot(2,1,2);imhist(X);title('histogram')
    Y = Bim_equalization(X);
    figure(2)
    subplot(2,1,1);imshow(Y); title('transformed image')
    subplot(2,1,2);imhist(Y);title('histogram')
end
