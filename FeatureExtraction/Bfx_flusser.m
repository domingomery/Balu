% [X,Xn] = Bfx_flusser(R,options)
% [X,Xn] = Bfx_flusser(R)
%
% Toolbox: Balu
%
%    Extract the four Flusser moments from binary image R.
%
%    options.show = 1 display mesagges.
%
%    X is a 4 elements vector:
%      X(i): Flusser-moment i for i=1,...,4.
%    Xn is the list of feature names.
%
%    Reference:
%    Sonka et al. (1998): Image Processing, Analysis, and Machine Vision,
%    PWS Publishing. Pacific Grove, Ca, 2nd Edition.
%
%    Example:
%      I = imread('testimg3.jpg');     % input image
%      R = Bim_segbalu(I);             % segmentation
%      [L,n] = bwlabel(R);             % regions
%      imshow(L,[])
%      X = [];
%      for i=1:n
%         [Xi,Xn] = Bfx_flusser(L==i); % Flusser moments
%         X = [X;Xi];
%      end
%      X
%
%   See also Bfx_standard, Bfx_hugeo, Bfx_fitellipse, Bfx_gupta.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [X,Xn] = Bfx_flusser(R,options)
if ~exist('options','var')
    options.show = 0;
end

if options.show == 1
    disp('--- extracting Flusser moments...');
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
u00 = (I0'*J0);
% u01 = (I0'*J1); not used
u02 = (I0'*J2);
u03 = (I0'*J3);
% u10 = (I1'*J0); not used
u20 = (I2'*J0);
u30 = (I3'*J0);
u11 = (I1'*J1);
u12 = (I1'*J2);
u21 = (I2'*J1);

II1 = (u20*u02-u11^2)/u00^4 ;
II2 = (u30^2*u03^2-6*u30*u21*u12*u03+4*u30*u12^3+4*u21^3*u03-3*u21^2*u12^2)/u00^10;
II3 = (u20*(u21*u03-u12^2)-u11*(u30*u03-u21*u12)+u02*(u30*u12-u21^2))/u00^7;
II4 = (u20^3*u03^2-6*u20^2*u11*u12*u03-6*u20^2*u02*u21*u03+9*u20^2*u02*u12^2 + 12*u20*u11^2*u21*u03+6*u20*u11*u02*u30*u03-18*u20*u11*u02*u21*u12-8*u11^3*u30*u03- 6*u20*u02^2*u30*u12+9*u20*u02^2*u21+12*u11^2*u02*u30*u12-6*u11*u02^2*u30*u21+u02^3*u30^2)/u00^11;

X  = [II1 II2 II3 II4];

Xn = [ 'Flusser-moment 1        '
       'Flusser-moment 2        '
       'Flusser-moment 3        '
       'Flusser-moment 4        '];

