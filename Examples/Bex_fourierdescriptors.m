% Bex_fourierdescriptors
%
% Toolbox: Balu
%    Example: Separation between T and Y
%
%    This example uses the eccentricity to separate T and Y.
%
%    The classification is performed using a threshold computed as the
%    average of the class centroids.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl

I = imread('letters.bmp');
J = imdilate(not(I),ones(3,3));
figure(1)
clf
imshow(J,[]);
L = bwlabel(J,4);
title('Select a charcter using mouse double-click...');
[j,i] = getpts;
R = L==L(fix(i(1)+0.5),fix(j(1)+0.5));

E = bwperim(R,4);

options.Nfourierdes = 60;
options.show = 0;

F = Bfx_fourierdes(R,options);

%[Ip,Jp] = find(E==1);               % Pixel of perimeter in (i,j)
% F= fourierdes(Ip,Jp);

figure(1);
imshow(E)

figure(3);
plot(log(abs(F)+1))

n = length(F);

Ff = zeros(1,n);
M = [2 3 4 5 6 7 8 10 12 14 16 18 20 22 25 30 35 40 50 60];
a = 0;
for i=1:20;
    subplot(4,5,i)
    m = M(i);
    Ff(1:m) = F(1:m);
    Ff(n-m+2:n) = F(n-m+2:n);
    f = ifft(Ff);
    f = [f f(1:2)];
    x = imag(f);
    y = -real(f);
    plot(x,y,'r')
    title(sprintf('m=%d',m))
    if (a==0)
        ax = axis;
        ax = [min(x)-20 max(x)+(2)+20 min(y)-20 max(y)+20];
        a = 1;
    end
    axis(ax);
    axis off
end
figure(4)
stem(abs(F))
title('Spectrum')