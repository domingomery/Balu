% [X,Xn,Xu] = Bfx_huint(I,R,options)
%
% Toolbox: Balu
%    Hu moments with intensity.
%
%    X is a 7 elements vector:
%      X(i): Hu-moment i for i=1,...,7.
%    Xn is the list of feature names (see Example to see how it works).
%
%    Reference: 
%    Hu, M-K.: "Visual Pattern Recognition by Moment Invariants",
%    IRE Trans. Info. Theory  IT-8:179-187: 1962.
%
%    Example:
%      I = imread('testimg1.jpg');           % input image
%      R = Bim_segbalu(I);                   % segmentation
%      J = double(I(:,:,1))/256;             % normalized red channel
%      options.show    = 1;                  % display results
%      [X,Xn] = Bfx_huint(J,R,options);      % Hu moments with intenisty
%      Bio_printfeatures(X,Xn)
%
%   See also Bfx_haralick, Bfx_clp, Bfx_fourier, Bfx_dct, Bfx_lbp.
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [X,Xn] = Bfx_huint(I,R,options)

if ~exist('options','var')
    options.show = 0;
end

if options.show == 1
    disp('--- extracting Hu moments with intensity...');
end
[Ireg,Jreg] = find(R==1);           % pixels in the region
im = mean(Ireg);
jm = mean(Jreg);
Kreg = R==1;
A    = length(Ireg);


I0  = ones(A,1);
J0  = ones(A,1);
I1 = Ireg - im*ones(A,1);
J1 = Jreg - jm*ones(A,1);
I2 = I1.*I1;
J2 = J1.*J1;
I3 = I2.*I1;
J3 = J2.*J1;
xreg = I(Kreg);
if (sum(xreg==0))
   xreg = ones(A,1);
end
J0X = double(J0).*double(xreg);

% Central moments
u00 = (I0'*J0X);
u002 = u00*u00;
u0025 = u00^2.5;
n02 = (I0'*J2)/u002;
n20 = (I2'*J0)/u002;
n11 = (I1'*J1)/u002;
n12 = (I1'*J2)/u0025;
n21 = (I2'*J1)/u0025;
n03 = (I0'*J3)/u0025;
n30 = (I3'*J0)/u0025;

f1 = n20+n02; 
f2 = (n20-n02)^2 + 4*n11^2;
f3 = (n30-3*n12)^2+(3*n21-n03)^2; 
f4 = (n30+n12)^2+(n21+n03)^2; 
f5 = (n30-3*n12)*(n30+n12)*((n30+n12)^2 - 3*(n21+n03)^2) + (3*n21-n03)*(n21+n03)*(3*(n30+n12)^2 - (n21+n03)^2);
f6 = (n20-n02)*((n30+n12)^2 - (n21+n03)^2) + 4*n11*(n30+n12)*(n21+n03);
f7 = (3*n21-n03)*(n30+n12)*((n30+n12)^2 - 3*(n21+n03)^2) - (n30-3*n12)*(n21+n03)*(3*(n30+n12)^2 - (n21+n03)^2);


X  = [f1 f2 f3 f4 f5 f6 f7];

Xn = [ 'Hu-moment-int 1         '
       'Hu-moment-int 2         '
       'Hu-moment-int 3         '
       'Hu-moment-int 4         '
       'Hu-moment-int 5         '
       'Hu-moment-int 6         '
       'Hu-moment-int 7         '];
