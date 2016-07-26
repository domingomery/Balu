% [X,Xn,] = Bfx_dct(I,R,options)
% [X,Xn] = Bfx_dct(I,options)
%
% Toolbox: Balu
%    DCT features
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
%      options.Ndct  = 64;                % imresize vertical
%      options.Mdct  = 64;                % imresize horizontal
%      options.mdct  = 2;                 % imresize frequency vertical
%      options.ndct  = 2;                 % imresize frequency horizontal
%      options.show    = 1;               % display results
%      I = imread('testimg1.jpg');        % input image
%      R = Bim_segbalu(I);                % segmentation
%      J = double(I(:,:,2))/256;          % normalized green channel
%      [X,Xn] = Bfx_dct(J,R,options);     % dct features
%      Bio_printfeatures(X,Xn)
%
%   See also Bfx_haralick, Bfx_clp, Bfx_gabor, Bfx_fourier, Bfx_lbp.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [X,Xn] = Bfx_dct(I,R,options)

if nargin==2;
    options = R;
    R = ones(size(I));
end

I(R==0) = 0;

N  = options.Ndct;
M  = options.Mdct;
n  = options.ndct;
m  = options.mdct;

N2 = round(N/2);
M2 = round(M/2);

if options.show
    disp('--- extracting dct features...');
end


Im = imresize(double(I),[N M]);
Fm = abs(dct2(Im));
F  = imresize(Fm(1:N2,1:M2),[n m]); 

LS = n*m;
X  = zeros(1,LS);
Xn = char(zeros(LS,24));
k = 0;
for i=1:n
    for j=1:m
        k = k + 1;
        s = sprintf('DCT(%d,%d)                  ',i,j);
        Xn(k,:)  = s(1:24);
        X(k) = F(i,j);
    end
end
