% [X,Xn] = Bfx_basicint(I,R,options)
% [X,Xn] = Bfx_basicint(I,options)
%
% Toolbox: Balu
%    Basic intensity features
%
%    X is the features vector, Xn is the list feature names (see Example to
%    see how it works).
%
%    Reference:
%    Kumar, A.; Pang, G.K.H. (2002): Defect detection in textured materials
%    using Gabor filters. IEEE Transactions on Industry Applications,
%    38(2):425-440.
%
%   Example:
%      options.mask  = 5;                      % Gauss mask for gradient computation
%      options.show  = 1;                      % display results
%      I = imread('testimg1.jpg');             % input image
%      R = Bim_segbalu(I);                     % segmentation
%      J = double(I(:,:,2))/256;               % normalized green channel
%      [X,Xn] = Bfx_basicint(J,R,options);     % basic intenisty features
%      Bio_printfeatures(X,Xn)
%
%   See also Bfx_haralick, Bfx_clp, Bfx_gabor, Bfx_fourier, Bfx_dct, Bfx_lbp.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [X,Xn] = Bfx_basicint(I,R,options)

if nargin==2;
    options = R;
    R = ones(size(I));
end

if ~isfield(options,'mask')
   options.mask = 15;
end

if options.show
    disp('--- extracting basic intensity features...');
end


E = bwperim(R,4);

ii = find(R==1);
jj = find(R==0, 1);
kk = E==1;

I = double(I);

I1 = Bim_d1(I,options.mask);
I2 = Bim_d2(I);

if ~isempty(jj)
   C = mean(abs(I1(kk)));
else
   C = -1;
end


J = I(ii);

G  = mean(J);
S  = std(J);
K  = kurtosis(J);
Sk = skewness(J);
D  = mean(I2(ii));

X = [G S K Sk D C];

Xn = [ 'Intensity Mean          '
       'Intensity StdDev        '
       'Intensity Kurtosis      '
       'Intensity Skewness      '
       'Mean Laplacian          '
       'Mean Boundary Gradient  '];

