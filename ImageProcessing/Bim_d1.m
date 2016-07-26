% function J = Bim_d1(I,m);
%
% Toolbox: Balu
%    First derivative of image X using a m x m Gauss operator.
%
%    Input data:
%       I grayvalue image.
%
%    Output:
%       J first derivative of I
%
%    Example: 
%       X = imread('testimg2.jpg');
%       I = rgb2gray(X);
%       figure(1)
%       imshow(I); title('original image')
%       J = Bimgd1(I,5);
%       figure(2)
%       imshow(J,[]); title('gradient')
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [Y,Yx,Yy] = Bim_d1(X,m)
sigma = m/8.5;
s2 = sigma^2;
Gx = zeros(m,m);
Gy = zeros(m,m);
c = (m-1)/2;
for i = 1:m
   x = i-c;
   x2 = (i-c)^2;
   for j=1:m
       y = j-c;
      y2 = (j-c)^2;
      ex = exp(-(x2+y2)/2/s2);
      Gx(i,j) = y*ex;
      Gy(i,j) = x*ex;
   end
end
mgx = sum(abs(Gx(:)))/2*(0.3192*m-0.3543);
Gx = Gx/mgx;
Gy = Gy/mgx;
Yx = conv2(X,Gx,'same');
Yy = conv2(X,Gy,'same');
Y0 = sqrt(Yx.*Yx+Yy.*Yy);
[N,M] = size(X);
Y = zeros(N,M);
Y(c:N-c,c:M-c) = Y0(c:N-c,c:M-c);
