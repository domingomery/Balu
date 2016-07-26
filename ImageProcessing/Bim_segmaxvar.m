% [R,E,J] = Bim_segmaxvar(I,p)
%
% Toolbox: Balu
%  Segmentation of an object with homogeneous background using as objective
%  function the maximization of the variance of a linear transformation of
%  the color image.
%
%  I: input image
%  p: threshold (default: p=-0.05) with p between -1 and 1.
%      A positive value is used to dilate the segmentation,
%      the negative to erode.
%  R: binary image of the object
%  E: binary image of the edge of the object
%  J: high contrast image of I.
%
%  Example:
%     I = imread('testimg1.jpg');
%     R = Bim_segmaxvar(I);
%     figure(1)
%     imshow(I); title('test image')
%     figure(2)
%     imshow(R); title('segmented image')
%
%     Repeat this examples for images testimg2, testimg3 and testimg4. Last
%     test image requires R = Bim_segbalu(I,-0.1) for better results.
%
% See also Bim_segmowgli, Bim_segotsu, Bim_segkmeans, Bio_segshow,
% Bim_segbalu.
%
% (c) G. Tejeda, PUC-DCC, May 2011

function [R,E,J] = Bim_segmaxvar(I,p)

I = double(I);
if size(I,3)==3
    k = Balgebraic(imresize(I,[64 64]));
    J = k(1)*I(:,:,1)+k(2)*I(:,:,2)+k(3)*I(:,:,3);
    J = Bim_maxmin(J);
else
    J = Bim_maxmin(I);
end
n = fix(size(J,1)/4);
if (mean2(J(1:n,1:n)) > 0.4)
    J = 1 - J;
end


t = graythresh(J);
if (~exist('p','var'))
    p = -0.05;
end
[R,E] = Bim_morphoreg(J,t+p);
end


function k=Balgebraic(RGB)

r     = RGB(:,:,1);
g     = RGB(:,:,2);
b     = RGB(:,:,3);
C     = cov([r(:) g(:) b(:)]);
vrr   = C(1,1);vgg=C(2,2);vbb=C(3,3);
vrg   = C(1,2);vrb=C(1,3);vgb=C(2,3);

w     = Broots([1 (vbb+vrr+vgg) (-vgb^2-vrg^2+vbb*vrr-vrb^2+vbb*vgg+vgg*vrr) -vrg^2*vbb+vbb*vgg*vrr+2*vgb*vrb*vrg-vrb^2*vgg-vgb^2*vrr]);

alpha = (vrb*vrg-vgb*vrr-vgb*w)./(vgb*vrg-vrb*vgg-vrb*w);
beta  = (-vrg^2+vgg*vrr+vgg*w+w*vrr+w.^2)./(vgb*vrg-vrb*vgg-vrb*w);

r1    = sqrt(1./(1+alpha.^2+beta.^2));
r2    = alpha.*r1;
r3    = beta.*r1;

J     = r1.^2*vrr+r2.^2*vgg+r3.^2*vbb+2*r1.*r2*vrg+2*r1.*r3*vrb+2*r2.*r3*vgb;

i     = J==max(J);

k=[r1(i) r2(i) r3(i)];

end

function y = Broots(x)

a    = x(1);
b    = x(2);
c    = x(3);
d    = x(4);

y    = zeros(1,3);
Q    = (2*b^3-9*a*b*c+27*a^2*d)^2-4*(b^2-3*a*c)^3;
if Q>0
    return;
end
Q    = sqrt(Q);
C    = (1/2*(Q+2*b^3-9*a*b*c+27*a^2*d))^(1/3);

y(1) = -b/(3*a)-C/(3*a)-(b^2-3*a*c)/(3*a*C);
y(2) = -b/(3*a)+C*(1+1i*sqrt(3))/(6*a)+(1-1i*sqrt(3))*(b^2-3*a*c)/(6*a*C);
y(3) = -b/(3*a)+C*(1-1i*sqrt(3))/(6*a)+(1+1i*sqrt(3))*(b^2-3*a*c)/(6*a*C);

y    = real(y);

end