% [X,Xn,options] = Bfx_bsif(I,R,options)
% [X,Xn,options] = Bfx_bsif(I,options)
% [X,Xn] = Bfx_bsif(I,R,options)
% [X,Xn] = Bfx_bsif(I,options)
%
% Toolbox: Balu
%    Binarized statistical image features
%
%    X is the features vector, Xn is the list of feature names (see Example
%    to see how it works).
%
%    It calculates the BSIF over the a regular grid of patches. The function
%    uses Juho Kannala and Esa Rahtu (see http://www.ee.oulu.fi/~jkannala/bsif/bsif.html).
%
%    It returns a matrix of uniform bsif descriptors for I, made by
%    concatenating histograms of each grid cell in the image.
%    Grid size is options.hdiv * options.vdiv
%
%    R is a binary image or empty. If R is given the bsif will be computed
%    the corresponding pixles R==0 in image I will be set to 0.
%
%     Output:                                            %%%%%%%%%%%%%%%%%%%%%revisar%%%%%%%%%%%%
%     X is a matrix of size ((hdiv*vdiv) x 59), each row has a
%         histogram corresponding to a grid cell. We use 59 bins.
%     options.x of size hdiv*vdiv is the x coordinates of center of ith grid cell
%     options.y of size hdiv*vdiv is the y coordinates of center of ith grid cell
%     Both coordinates are calculated as if image was a square of side length 1.
%
%     References:
%     J. Kannala and E. Rahtu. Bsif: Binarized statistical image features.  
%     In Pattern Recognition (ICPR), 2012 21st International Conference on,
%     pages 1363--1366. IEEE, 2012.
%
%
%     Example 1:
%      options.vdiv = 1;                  % one vertical divition
%      options.hdiv = 1;                  % one horizontal divition
%      options.filter  = 7;               % use filter 7x7
%      options.bits  = 11;                % use 11 bits filter
%      options.mode  = 'h';               % return histogram
%      I = imread('testimg1.jpg');        % input image
%      J = I(120:219,120:239,2);          % region of interest (green)
%      figure(1);imshow(J,[])             % image to be analyzed
%      [X,Xn] = Bfx_bsif(J,[],options);   % BSIF features
%      figure(2);bar(X)                   % histogram
%
%     Example 2:
%      options.vdiv = 1;                  % one vertical divition
%      options.hdiv = 1;                  % one horizontal divition
%      options.filter  = 7;               % use filter 7x7
%      options.bits  = 11;                % use 11 bits filter
%      options.mode  = 'nh';              % return normilized histogram
%      I = imread('testimg1.jpg');        % input image
%      J = I(120:219,120:239,2);          % region of interest (green)
%      figure(1);imshow(J,[])             % image to be analyzed
%      [X,Xn] = Bfx_bsif(J,[],options);   % BSIF features
%      figure(2);bar(X)                   % histogram
%
%     Example 3:
%      options.vdiv = 1;                  % one vertical divition
%      options.hdiv = 1;                  % one horizontal divition
%      options.filter  = 7;               % use filter 7x7
%      options.bits  = 11;                % use 11 bits filter
%      options.mode  = 'im';              % return image represetation
%      %(image is only aviable for vdiv = 1 y hdiv = 1)
%      I = imread('testimg1.jpg');        % input image
%      J = I(120:219,120:239,2);          % region of interest (green)
%      figure(1);imshow(J,[])             % image to be analyzed
%      [X,Xn] = Bfx_bsif(J,[],options);   % BSIF features
%      figure(2);imshow(X,[]);            % display image
%
%
%   See also Bfx_lbp, Bfx_gabor, Bfx_clp, Bfx_fourier, Bfx_dct.
%
% (c) Erick Svec
%
function [X,Xn] = Bfx_bsif(I,R,options)

if nargin==2;
    options = R;
    R = ones(size(I));
end

vdiv = options.vdiv;
hdiv = options.hdiv;

if ~isfield(options,'show')
    options.show = 0;
end

if options.show == 1
    disp('--- extracting binarized statistical image features...');
end


filename=['ICAtextureFilters_' int2str(options.filter) 'x' int2str(options.filter) '_' int2str(options.bits) 'bit'];
load(filename, 'ICAtextureFilters');


if ~isempty(R);
    I(R==0) = 0;
end

if strcmp(options.mode,'im')
    if vdiv == 1 && hdiv == 1
        X = bsif(I,ICAtextureFilters,options.mode);
        Xn = options;
        return
    else
        throw(MException('MATLAB:odearguments:InconsistentDataType', 'Invalid Options Set: vdiv and hdiv must be 1 for mode im (when try retrieve image representation)'));
        return
    end
end
[N,M] = size(I);
w = N/vdiv;
h = M/hdiv;

X = [];
Xn = options;

for i = 1:vdiv
    for j = 1:hdiv
        part = I(i*w-w+1:w*i,j*h-h+1:h*j);
        code_img = bsif(part,ICAtextureFilters,options.mode);
        X = [X code_img];
    end
end

end