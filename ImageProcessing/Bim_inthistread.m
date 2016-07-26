% function h = Bim_inthistread(H,i1,j1,i2,j2,options)
%
%
% Toolbox: Balu
%    Histogram of a part of an image using integral histograms.
%
%    Input data:
%       H integral histogram.
%       (i1,j1,i2,j2) rectangle of the image
%
%    Output:
%       h histogram
%
%    Example:
%       I = imread('rice.png');
%       H = Bim_inthist(I,256);
%       i1=70;j1=50;i2=130;j2=190;
%       K = I(i1:i2,j1:j2);
%       t = hist(K(:),1:256);
%       h = Bim_inthistread(H,i1,j1,i2,j2);
%       compare(h,t)
%
%    See also Bim_inthist.
%
% (c) D.Mery, PUC-DCC, 2012
% http://dmery.ing.puc.cl
%

function h = Bim_inthistread(H,x1,y1,x2,y2,n)
B = size(H,3);

if ~exist('n','var')
    n = 1;
end

if n==1
    h = suminthist(H,x1,y1,x2,y2,B);
else
    %[x1 y1 x2 y2]
    N  = size(H,1);
    M  = size(H,2);
    t = 0;
    h = zeros(1,n^2*B);
    dx = (x2+1-x1)/n;
    dy = (y2+1-y1)/n;
    xx1 = x1;
    for i=1:n
        xx2 = min([round(xx1+dx)-1 N]);
        yy1 = y1;
        for j=1:n
            t = t+1;
            yy2 = min([round(yy1+dy)-1 M]);
            % [i j xx1 yy1 xx2 yy2]
            h(:,indices(t,B)) = suminthist(H,xx1,yy1,xx2,yy2,B);
            yy1 = yy2+1;
        end
        xx1 = xx2+1;
    end

% disp('****')

end
end
function h = suminthist(H,x1,y1,x2,y2,B)
h = zeros(1,B);
t = zeros(1,B);
% [0 0 x1 y1 x2 y2]
% size(H)

h(:) = H(x2,y2,:);


if and((x1>1),(y1>1))
    t(:) = H(x1-1,y1-1,:);
    h(:) = h(:)+t(:);
end
if x1>1
    t(:) = H(x1-1,y2,:);
    h(:) = h(:)-t(:);
end
if y1>1
    t(:) = H(x2,y1-1,:);
    h(:) = h(:)-t(:);
end

end


