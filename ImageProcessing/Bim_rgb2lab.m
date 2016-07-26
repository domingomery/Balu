% K = Bim_rgb2lab(I,M)
%
% Toolbox: Balu
%    Conversion RGB->L*a*b according to parameters M estimated with 
%    function Bim_labparam.m after
%    Leon,K.; Mery,D.; Pedreschi,F.;Leon,J.(2006): Color measurement in
%    L*a*b* units from RGB digital images. Food Research International,
%    39(10):1084-1091.
%  
%    I is RGB image
%    M parameters of the conversion
%    K is L*a*b* image
%
% (c) D.Mery, PUC-DCC, May. 2008
% http://dmery.ing.puc.cl

function K = Bim_rgb2lab(I,M)

r = double(I(:,:,1));
g = double(I(:,:,2));
b = double(I(:,:,3));

R = r(:);
G = g(:);
B = b(:);
n = length(R);

if (size(M,2)==4) % linear model
    X = [R G B ones(n,1)];
else              % quadratic model
    X = [R G B R.*G R.*B G.*B R.*R G.*G B.*B ones(n,1)];
end
L = X*M(1,:)';
A = X*M(2,:)';
B = X*M(3,:)';
   
l = zeros(size(r));
a = l;
b = l;

l(:) = L;
a(:) = A;
b(:) = B;

K = zeros(size(I));

K(:,:,1) = Bim_sat(l,0,100);
K(:,:,2) = Bim_sat(a,-120,120);
K(:,:,3) = Bim_sat(b,-120,120);
