% [X,Xn] = Bfx_hog(I,options)
%
% Toolbox: Balu
%    Histogram of Orientated Gradients features
%
%    X is the features vector, Xn is the list of feature names (see Example
%    to see how it works).
%
%    options.nj;  : number of HOG windows per bound box
%    options.ni   : in i (vertical) and j (horizaontal) direction
%    options.B    : number of histogram bins
%    options.show : show histograms (glyphs)
%
%     Example:
%        options.nj    = 20;             % 10 x 20 
%        options.ni    = 10;             % histograms
%        options.B     = 9;              % 9 bins
%        options.show  = 1;              % number of neighbor samples
%        I = imread('testimg1.jpg');     % input image
%        J = rgb2gray(I);
%        figure(1);imshow(J,[]);
%        figure(2);
%        [X,Xn] = Bfx_hog(J,options);    % HOG features (see gradients
%                                        % arround perimeter.
%
%   See also Bfx_phog, Bfx_lbp.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl
%

function [X,Xn,options] = Bfx_lbphogi(I,R,options)

if nargin==2;
    options = R;
    R = [];
end

[X_lbp,Xn_lbp] = Bfx_lbpi(I,R,options);
[X_hog,Xn_hog] = Bfx_hogi(I,R,options);

X  = [X_lbp  X_hog ];
Xn = [Xn_lbp; Xn_hog];

