% [X,Xn] = Bfx_hog(I,options)
%
% Toolbox: Balu
%    Histogram of Orientated Gradients features
%
%    X is the features vector, Xn is the list of feature names (see Example
%    to see how it works).
%
%    options.nj;  : number of HOG windows per bound box
%    options.ni   : in i (vertical) and j (horizaontal) direction
%    options.B    : number of histogram bins
%    options.show : show histograms (glyphs)
%
%     Example:
%        options.nj    = 20;             % 10 x 20 
%        options.ni    = 10;             % histograms
%        options.B     = 9;              % 9 bins
%        options.show  = 1;              % show results
%        I = imread('testimg1.jpg');     % input image
%        J = rgb2gray(I);
%        figure(1);imshow(J,[]);
%        figure(2);
%        [X,Xn] = Bfx_hog(J,options);    % HOG features (see gradients
%                                        % arround perimeter).
%
%   See also Bfx_phog, Bfx_lbp.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl
%

function [X,Xn,options] = Bfx_hog(I,R,options)

if nargin==2;
    options = R;
    R = ones(size(I));
end

if ~isfield(options,'normalize')
    options.normalize = 0;
end


nj    = options.nj;       % number of HOG windows per bound box
ni    = options.ni;       % in i (vertical) and j (horizaontal) direction
B     = options.B;        % number of histogram bins
show  = options.show;     % show histograms

N     = size(I,1);
M     = size(I,2);
X     = zeros(1,nj*ni*B); % column vector with zeros
I     = double(I);
dj    = floor(M/(nj+1));
di    = floor(N/(ni+1));
t     = 0;
hj    = [-1,0,1];
hi    = -hj';
Gj    = imfilter(I,hj);
Gi    = imfilter(I,hi);
A     = atan2(Gi,Gj);
A     = mod(A,pi);
G     = ((Gi.^2)+(Gj.^2)).^.5;
K     = 4*di*dj;
ang   = zeros(K,1);
mag   = zeros(K,1);
if show
    J     = zeros(N,M);
    ss    = min([N/ni M/nj])*0.40;
end
w = zeros(ni,nj,B);
for i = 1:ni
    ii = (i-1)*di+1:(i+1)*di;
    i0 = mean(ii);
    for j = 1:nj
        jj     = (j-1)*dj+1:(j+1)*dj;
        j0     = mean(jj);
        t      = t+1;
        ang(:) = A(ii,jj);
        mag(:) = G(ii,jj);
        X2     = zeros(B,1);
        for b = 1:B
            q      = find(ang<=pi*b/B);
            X2(b)  = X2(b)+sum(mag(q));
            ang(q) = 5;
        end
        X2                = X2/(norm(X2)+0.01);
        X(1,indices(t,B)) = X2;
        % w(i,j,:)          = X2;
        if show
            for b=1:B
                alpha = pi*b/B;                
                q     = -ss:ss;
                qi    = round(i0+q*cos(alpha));
                qj    = round(j0+q*sin(alpha));
                qk    = qi+(qj-1)*N;
                J(qk) = J(qk)+X2(b);
            end
            J(round(i0),round(j0))=1;
        end
    end
end
if show
    imshow(round(J/max2(J)*256),jet)
    options.J = J;
end
options.w = w;
Xn = char(zeros(nj*ni*B,24));
Xn(:,1) = 'H'; Xn(:,2) = 'O'; Xn(:,3) = 'G';

J = uint8(zeros(N,M,3));
G(G>363) = 363; % max2(Gi) = max2(Gi) = 256 => max2(G) = sqrt(2*256^2)
J(:,:,1) = uint8(round(G/363*255));
J(:,:,2) = uint8(round(A/pi*255));
options.Ihog = J;
if options.normalize
    X = X/sum(X);
end




