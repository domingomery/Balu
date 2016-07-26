% [X,Xn] = Bfx_vlhog(I,options)
%
% Toolbox: Balu
%    Histogram of Orientated Gradients features using Vlfeat Toolbox.
%
%    X is the features vector, Xn is the list of feature names (see Example
%    to see how it works).
%
%    options.cellsize  : size of the cells in pixels
%    options.variant   : 1 for UoCTTI, 2 for Dalal-Triggs
%    options.show      : 1 shows oriented histograms
%
%     Example:
%        options.cellsize = 32;          % 32 x 32 
%        options.variant  = 1;           % UoCTTI
%        options.show  = 1;              % show results
%        I = imread('testimg1.jpg');     % input image
%        J = rgb2gray(I);
%        figure(1);imshow(J,[]);
%        figure(2);
%        [X,Xn] = Bfx_vlhog(J,options);  % HOG features (see gradients
%                                        % arround perimeter).
%
%   See also Bfx_phog, Bfx_lbp, Bfx_hog, vl_hog.
%
% (c) GRIMA-DCCUC, 2012
% http://grima.ing.puc.cl
%

function [X,Xn,options] = Bfx_vlhog(I,R,options)

if nargin==2;
    options = R;
    R = ones(size(I));
end

if ~isfield(options,'variant')
    options.varian = 1;
end

if ~isfield(options,'show')
    options.show = 1;
end

if options.variant == 1
    varname = 'UoCTTI';
else
    varname = 'DalalTriggs';
end

options.hog = vl_hog(im2single(I),options.cellsize,'variant',varname);
if options.show==1
    figure
    options.Ir = vl_hog('render',options.hog);
    imshow(options.Ir,[]);
end
X = options.hog(:)';
n = length(X);
Xn = zeros(n,24);







