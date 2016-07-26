% [X,Xn] = Bfx_gabor(I,R,options)
% [X,Xn] = Bfx_gabor(I,options)
%
% Toolbox: Balu
%    Gabor features
%
%    X is the features vector, Xn is the list of feature names (see Example 
%    to see how it works).
%
%    Reference:
%    Kumar, A.; Pang, G.K.H. (2002): Defect detection in textured materials
%    using Gabor filters. IEEE Transactions on Industry Applications,
%    38(2):425-440.
%
%   Example:
%      options.Lgabor  = 8;                 % number of rotations
%      options.Sgabor  = 8;                 % number of dilations (scale)
%      options.fhgabor = 2;                 % highest frequency of interest
%      options.flgabor = 0.1;               % lowest frequency of interest
%      options.Mgabor  = 21;                % mask size
%      options.show    = 1;                 % display results
%      I = imread('testimg1.jpg');          % input image
%      R = Bim_segbalu(I);                  % segmentation
%      J = I(:,:,2);                        % green channel
%      [X,Xn] = Bfx_gabor(J,R,options);  % Gabor features
%      Bio_printfeatures(X,Xn)
%
%   See also Bfx_haralick, Bfx_clp, Bfx_fourier, Bfx_dct, Bfx_lbp.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [X,Xn] = Bfx_gabor(I,R,options)

if nargin==2;
    options = R;
    R = ones(size(I));
end

L  = options.Lgabor;
S  = options.Sgabor;
fh = options.fhgabor;
fl = options.flgabor;
M  = options.Mgabor;

if options.show
    disp('--- extracting Gabor features...');
end

alpha = (fh/fl)^(1/(S-1));
sx = sqrt(2*log(2))*(alpha+1)/2/pi/fh/(alpha-1);
sy = sqrt(2*log(2)-(2*log(2)/2/pi/sx/fh)^2)/(2*pi*tan(pi/2/L)*(fh-2*log(1/4/pi^2/sx^2/fh)));
u0 = fh;

k = R==1;

g = zeros(S,L);
size_out = size(I)+[M M]-1;
Iw = fft2(I,size_out(1),size_out(2));
n1 = (M+1)/2;
[NN,MM] = size(I);

for p=1:S;
    for q=1:L
        f = Bgabor_pq(p,q,L,sx,sy,u0,alpha,M);
        %This convolution is very slow:
        %Ir = conv2(I,real(f),'same');
        %Ii = conv2(I,imag(f),'same');
        %Iout = sqrt(Ir.*Ir + Ii.*Ii);
        Ir = real(ifft2(Iw.*fft2(real(f),size_out(1),size_out(2))));
        Ii = real(ifft2(Iw.*fft2(imag(f),size_out(1),size_out(2))));
        Ir = Ir(n1:n1+NN-1,n1:n1+MM-1);
        Ii = Ii(n1:n1+NN-1,n1:n1+MM-1);
        Iout = sqrt(Ir.*Ir + Ii.*Ii);
        g(p,q) = mean(Iout(k));
    end;
end
gmax = max(g(:));
gmin = min(g(:));
J = (gmax-gmin)/gmin;
X = [g(:); gmax; gmin; J]';

LS = L*S;
Xn = char(zeros(LS+3,24));
k = 0;
for p=1:S;
    for q=1:L
        k = k + 1;
        Xn(k,:) = sprintf('Gabor(%d,%d)              ',p,q);
    end
end

Xn(LS+1,:)  = 'Gabor-max               ';
Xn(LS+2,:)  = 'Gabor-min               ';
Xn(LS+3,:)  = 'Gabor-J                 ';

% f = Bgabor_pq(p,q,L,S,sx,sy,u0,alpha,M)
%
% Toolbox: Balu
%    Gabor kernel. See details in:
%    Kumar, A.; Pang, G.K.H. (2002): Defect detection in textured materials
%    using Gabor filters. IEEE Transactions on Industry Applications,
%    38(2):425-440.
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl
%

function f = Bgabor_pq(p,q,L,sx,sy,u0,alpha,M)

f = zeros(M,M);
sx2 = sx*sx;
sy2 = sy*sy;
c = (M+1)/2;
ap = alpha^-p;
tq = pi*(q-1)/L;

f_exp = 2*pi*sqrt(-1)*u0;
for i=1:M
    x = i - c;
    for j=1:M
        y = j - c;
        x1 = ap*(x*cos(tq)+y*sin(tq));
        y1 = ap*(y*cos(tq)-x*sin(tq));
        f(i,j) = exp(-0.5*(x1*x1/sx2+y1*y1/sy2))*exp(f_exp*x1);
    end
end
f = ap*f/2/pi/sx/sy;