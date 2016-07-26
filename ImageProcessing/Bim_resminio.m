% Fs = Bim_resminio(G,h,method)
%
% Toolbox: Balu
%  Image restoration using MINIO criterium.
%
%  G : blurred image
%  h : PSF
%  method = 1: minio ||f_N-g|| -> min (default)
%         = 2: ||f|| -> min
%  Fs: restored image
%
%  See details in:
%  Mery, D.; Filbert, D. (2000): A Fast Non-iterative Algorithm for the
%  Removal of Blur Caused by Uniform Linear Motion in X-ray Images. In
%  Proceedings of the 15th World Conference on Non-Destructive Testing,
%  Oct. 15-21, Roma.
%
%  Example:
%     F = double(imread('circuit.tif'));
%     n = 45;
%     h = ones(1,n)/n;
%     G = conv2(F,h,'valid');
%     Fs = Bim_resminio(G,h);
%     figure(1)
%     imshow(F,[]); title('original image')
%     figure(2)
%     imshow(G,[]); title('blurred image')
%     figure(3)
%     imshow(Fs,[]); title('restored image')
%
% D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function Fs = Bim_resminio(G,h,method)
if ~exist('method','var')
    method = 1;
end
n = length(h);
N = size(G,2);
M = N + n -1;
H = zeros(N,M);
for i=1:N
    H(i,i:i+n-1) = h;
end
lambda = 1e6;
if method == 1
    P  = [eye(N,N) zeros(N,n-1)];
    Fs = G*(lambda*H+P)*(inv(lambda*H'*H+P'*P))';
else
    Fs = G*H*(inv(lambda*H'*H+eye(M,M)))';
end

