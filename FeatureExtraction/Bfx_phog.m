% [X,Xn,Xu] = Bfx_phog(I,R,options)
% [X,Xn,Xu] = Bfx_phog(I,options)
%
% Toolbox: Balu
%
%    Pyramid Histogram of Oriented Gradients based on implementation by
%    Anna Bosch from
%
%    http://www.robots.ox.ac.uk/~vgg/research/caltech/phog.html
%
%   IN:
%	I   - Images of size MxN (Color or Gray)
%	options.bin - Number of bins on the histogram
%   options.L - number of pyramid levels
%
%   OUT:
%    X is the features vector, Xn is the list of feature names (see Example 
%    to see how it works).
%
%   Reference:
%   Dalal, N. & Triggs, B. (2005): Histograms of oriented gradients for
%   human detection, Proceedings of the Conference on Computer Vision and 
%   Pattern Recognition, Vol. 1, 886-893
%
%   Example:
%      options.bin  = 9;                  % bins on the histogram
%      options.L    = 3;                  % pyramides levels
%      options.show    = 1;               % display results
%      I = imread('testimg1.jpg');        % input image
%      J = double(I(:,:,2))/256;          % normalized green channel
%      [X,Xn] = Bfx_phog(J,options);      % phog features
%      Bio_printfeatures(X,Xn)
%
%   See also Bfx_haralick, Bfx_clp, Bfx_gabor, Bfx_fourier, Bfx_lbp.
%
% D.Mery, A. Soto PUC-DCC, Jun. 2010
% http://dmery.ing.puc.cl


function [X,Xn] = Bfx_phog(I,R,options)

if nargin==2;
    options = R;
    R = ones(size(I));
end

I(R==0) = 0;

bin = options.bin;
L   = options.L;

if options.show
    disp('--- extracting phog features...');
end

roi   = [1 size(I,1) 1 size(I,2)]';
angle = 360;

if size(I,3) == 3
    G = rgb2gray(I);
else
    G = I;
end

if sum(sum(G))>100
    E = edge(G,'canny');
    [GradientX,GradientY] = gradient(double(G));
    Gr = sqrt((GradientX.*GradientX)+(GradientY.*GradientY));
    
    index = GradientX == 0;
    GradientX(index) = 1e-5;
    
    A = ((atan2(GradientY,GradientX)+pi)*180)/pi;
    
    [bh bv] = BphogbinMatrix(A,E,Gr,angle,bin);
else
    bh = zeros(size(I,1),size(I,2));
    bv = zeros(size(I,1),size(I,2));
end

bh_roi = bh(roi(1,1):roi(2,1),roi(3,1):roi(4,1));
bv_roi = bv(roi(1,1):roi(2,1),roi(3,1):roi(4,1));
X = BphogDescriptor(bh_roi,bv_roi,L,bin)';
n = length(X);

Xn = char(zeros(n,24));
for i=1:n
    s = sprintf('phog(%d)                   ',i);
    Xn(i,:)  = s(1:24);
end

end



function p = BphogDescriptor(bh,bv,L,bin)
% anna_PHOGDESCRIPTOR Computes Pyramid Histogram of Oriented Gradient over a ROI.
%
%    Pyramid Histogram of Oriented Gradients based on implementation by
%    Anna Bosch from
%
%    http://www.robots.ox.ac.uk/~vgg/research/caltech/phog.html
%
%IN:
%	bh - matrix of bin histogram values
%	bv - matrix of gradient values
%   L - number of pyramid levels
%   bin - number of bins
%
%OUT:
%	p - pyramid histogram of oriented gradients (phog descriptor)

p = [];
for b=1:bin
    ind = bh==b;
    p = [p;sum(bv(ind))];
end

cella = 1;
for l=1:L
    x = fix(size(bh,2)/(2^l));
    y = fix(size(bh,1)/(2^l));
    xx=0;
    yy=0;
    while xx+x<=size(bh,2)
        while yy +y <=size(bh,1)
            
            bh_cella = bh(yy+1:yy+y,xx+1:xx+x);
            bv_cella = bv(yy+1:yy+y,xx+1:xx+x);
            
            for b=1:bin
                ind = bh_cella==b;
                p = [p;sum(bv_cella(ind))];
            end
            yy = yy+y;
        end
        cella = cella+1;
        yy = 0;
        xx = xx+x;
    end
end
if sum(p)~=0
    p = p/sum(p);
end
end


function [bm bv] = BphogbinMatrix(A,E,G,angle,bin)
% anna_BINMATRIX Computes a Matrix (bm) with the same size of the image where
% (i,j) position contains the histogram value for the pixel at position (i,j)
% and another matrix (bv) where the position (i,j) contains the gradient
% value for the pixel at position (i,j)
%
%    Pyramid Histogram of Oriented Gradients based on implementation by
%    Anna Bosch from
%
%    http://www.robots.ox.ac.uk/~vgg/research/caltech/phog.html
%
%                
%IN:
%	A - Matrix containing the angle values
%	E - Edge Image
%   G - Matrix containing the gradient values
%	angle - 180 or 360%   
%   bin - Number of bins on the histogram 
%	angle - 180 or 360
%OUT:
%	bm - matrix with the histogram values
%   bv - matrix with the graident values (only for the pixels belonging to
%   and edge)

[contorns,n] = bwlabel(E);  
X = size(E,2);
Y = size(E,1);
bm = zeros(Y,X);
bv = zeros(Y,X);

nAngle = angle/bin;

for i=1:n
    [posY,posX] = find(contorns==i);    
    for j=1:size(posY,1)
        pos_x = posX(j,1);
        pos_y = posY(j,1);
        
        b = ceil(A(pos_y,pos_x)/nAngle);
        if b==0, bin= 1; end
        if G(pos_y,pos_x)>0
            bm(pos_y,pos_x) = b;
            bv(pos_y,pos_x) = G(pos_y,pos_x);                
        end
    end
end
end

