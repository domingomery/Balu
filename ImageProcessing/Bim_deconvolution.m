% Fs = Bim_deconvolution(G,h,C,SNR2,a)
%
% Toolbox: Balu
%  Image restoration using inverse filtering.
%
%                           (  H*(u,v) )a (       H*(u,v)')     )1-a
% Inverse Filter = W(u,v) = (----------)  (---------------------)
%                           (|H(u,v)|^2)  (|H(u,v)|^2 + C/SNR^2 )
%
% SNR^2 = 1/(r^rho)                  r = sqrt(u^2+v^2)
%
%  G : blurred image
%  h : PSF
%
%  See details in:
%  Castleman, K (1996): Digital Image Processing, Pretince Hall.
%
%  Example:
%  I = imread('saturn.png'); F = double(rgb2gray(imresize(I,[300 240])));
%  figure(1);imshow(F,[]); title('original image F')
%  n = 11; h = ones(n,n)/n/n;
%  G = conv2(F,h,'valid'); Gr = G+2*randn(size(G));
%  figure(2);imshow(Gr,[]); title('degraded image with noise')
%  Fs = Bim_deconvolution(G,h,0.025,0.01);
%  figure(3);imshow(Fs,[]); title('restored image')
%
%
% D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [Fs,SNR2] = Bim_deconvolution(G,h,C,SNR2,a)


if not(exist('C','var'))
    C=0;
end

if not(exist('a','var'))
    a=0;
end


[N,M]               = size(G);
[n,m]               = size(h);

N2                  = 2*(N+n-1);
M2                  = 2*(M+m-1);


Ge                  = zeros(N2,M2);
Ge(n:N+n-1,m:M+m-1) = G;
He                  = zeros(N2,M2);
He(1:n,1:m)         = h;
FG                  = fft2(Ge);
FH                  = fft2(He);
AFH                 = abs(FH);
AFH2                = AFH.*AFH;

if not(exist('SNR2','var'))
    SNR2 = ones(N2,M2);
else

    if length(SNR2(:)==1) %#ok<ISMT>
        rho = SNR2;
        SNR2 = zeros(N2,M2);
        % SNR2(r) = 1/(r^rho)
        for i=1:N2
            for j=1:M2
                SNR2(i,j) = 1/(sqrt((i-N2/2)^2+(j-M2/2)^2)^rho);
            end
        end
        %imshow(SNR2(1:2:end,1:2:end),[])
        %title('SNR^2')
        %disp('presione enter...')
        %pause
    end
end




FFa = (conj(FH)./AFH2).^a;

FFs = FFa.*(conj(FH)./(AFH2 + C./SNR2).^(1-a)).*FG;

Fs = abs(ifft2(FFs));

Ns = N+n-1;
Ms = M+m-1;
k=1;
Fs = Fs(k:k+Ns-1,k:k+Ms-1);
