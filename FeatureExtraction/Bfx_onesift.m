% [X,Xn,options] = Bfx_onesift(I,R,options)
% [X,Xn,options] = Bfx_onesift(I,options)
% [X,Xn] = Bfx_onesift(I,R,options)
% [X,Xn] = Bfx_onesift(I,options)
%
% Toolbox: Balu
%    Extract only one SIFT descriptor of region R of image I.
%
%    X is the features vector, Xn is the list of feature names (see Example
%    to see how it works).
%
%    R is a binary image or empty. If R is given the sift will be computed
%    in the region defined by the piexels where R==1.
%
%     Output:
%
%     References:
%     D. G. Lowe, Distinctive image features from scale-invariant
%        keypoints. IJCV, vol. 2, no. 60, pp. 91-110, 2004.
%
%     A. Vedaldi, B. Fulkerson: VLFeat: An Open and Portable Library
%        of Computer Vision Algorithms, 2008 (http://www.vlfeat.org/)
%
%
%     Example 5:
%      options.vdiv = 1;                  % one vertical divition
%      options.hdiv = 1;                  % one horizontal divition
%      options.semantic = 1;              % semantic LBP
%      options.samples = 8;               % number of neighbor samples
%      options.sk      = 0.25;            % angle sampling
%      options.weight  = 9;               % angle sampling
%      I = imread('testimg1.jpg');        % input image
%      J = I(120:219,120:239,2);          % region of interest (green)
%      [X,Xn] = Bfx_lbp(J,[],options);    % weighted LBP features
%      bar(X)                             % histogram
%   See also Bfx_gabor, Bfx_clp, Bfx_fourier, Bfx_dct.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl
%
function [X,Xn,options] = Bfx_onesift(I,R,options)

%if nargin==2;
%    options = R;
%    R = ones(size(I));
%end

if isempty(R)
    R = ones(size(I));
end

[ii,jj] = find(R==1);
ff = [mean(jj); mean(ii); sqrt(length(ii))/2; 0];

[ff,dd] = vl_sift(single(I),'Frames',ff);

X = dd';

n = 128;
Xn = char(zeros(n,24));
for k=1:n
        s = sprintf('SIFT(%d)                            ',k);
        Xn(k,:) = s(1:24);
end



