% [X,Xn] = Bfx_hugeo(R)
% [X,Xn] = Bfx_hugeo(R,options)
%
% Toolbox: Balu
%
%    Extract the seven Hu moments from binary image R.
%
%    options.show = 1 display mesagges.
%
%    X is a 7 elements vector:
%      X(i): Hu-moment i for i=1,...,7.
%      Xn is the list of feature names.
%
%    Reference: 
%    Hu, M-K.: "Visual Pattern Recognition by Moment Invariants",
%    IRE Trans. Info. Theory  IT-8:179-187: 1962.
%
%    Example:
%      I = imread('testimg3.jpg');     % input image
%      R = Bim_segbalu(I);             % segmentation
%      [L,n] = bwlabel(R);             % regions
%      imshow(L,[])
%      X = [];
%      for i=1:n
%         [Xi,Xn] = Bfx_hugeo(L==i);   % Hu moments
%         X = [X;Xi];
%      end
%      X
%
%   See also Bfx_basicgeo, Bfx_gupta, Bfx_fitellipse, Bfx_flusser.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [X,Xn] = Bfx_hugeo(R,options)
if ~exist('options','var')
    options.show = 0;
end

if options.show == 1
    disp('--- extracting Hu moments...');
end
[Ireg,Jreg] = find(R==1);           % pixels in the region

i_m = mean(Ireg);
j_m = mean(Jreg);
A   = length(Ireg);
I0  = ones(A,1);
J0  = ones(A,1);
I1 = Ireg - i_m*ones(A,1);
J1 = Jreg - j_m*ones(A,1);
I2 = I1.*I1;
J2 = J1.*J1;
I3 = I2.*I1;
J3 = J2.*J1;
% Central moments
u00 = A; %u00 = m00 = (I0'*J0);
u002 = u00*u00;
u0025 = u00^2.5;
% u0015 = u00^1.5; not used
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

Xn = [ 'Hu-moment 1             '
       'Hu-moment 2             '
       'Hu-moment 3             '
       'Hu-moment 4             '
       'Hu-moment 5             '
       'Hu-moment 6             '
       'Hu-moment 7             '];

