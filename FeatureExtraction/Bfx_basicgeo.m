% [X,Xn] = Bfx_geobasic(R,options)
%
% Toolbox: Balu
%
%    Standard geometric features of a binary image R. This function calls
%    regionprops of Image Processing Toolbox.
%
%    options.show = 1 display mesagges.
%
%    X is the feature vector
%    Xn is the list of feature names.
%
%    See expamle:
%
%   Example:
%      I = imread('testimg1.jpg');            % input image
%      R = Bim_segbalu(I);                    % segmentation
%      [X,Xn] = Bfx_basicgeo(R);              % basic geometric features
%      Bio_printfeatures(X,Xn)
%
%   See also Bfx_fitellipse, Bfx_hugeo, Bfx_gupta, Bfx_flusser.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [X,Xn] = Bfx_basicgeo(R,options)

if ~exist('options','var')
    options.show = 0;
end



N = size(R,1);


warning off
stats = regionprops(uint8(R),'EulerNumber','ConvexArea','EquivDiameter','Solidity','MajorAxisLength','Extent','MinorAxisLength','Orientation','FilledArea','Eccentricity');
warning on
E = bwperim(R,4);

% center of gravity
[Ireg,Jreg] = find(R==1);           % pixels in the region
Kreg        = Ireg+Jreg*N-N;        % pixels of region stored in a vector
i_m         = mean(Ireg);           % center of gravity in i direction
j_m         = mean(Jreg);           % center of gravity in i direction

% standard features
% Perimeter
if options.show == 1
    disp('--- extracting standard geometric features...');
end
L8 = sum(sum(bwperim(R,8)));
L4 = sum(sum(bwperim(R,4)));
L = (3*L4+L8)/4;

% Area
A = bwarea(R);

% Roundness
Roundness = 4*A*pi/L^2;

% heigh & width
i_max       = max(Ireg);
i_min       = min(Ireg);
j_max       = max(Jreg);
j_min       = min(Jreg);
height      = i_max-i_min+1;        % height
width       = j_max-j_min+1;        % width

% Danielsson shape factor (see Danielsson, 1977)
TD = double(bwdist(not(R),'chessboard'));
dm = mean(TD(Kreg));
Gd          = A/9/pi/dm^2;

X = [
    i_m
    j_m
    height
    width
    A
    L
    Roundness
    Gd
    stats.EulerNumber
    stats.EquivDiameter
    stats.MajorAxisLength
    stats.MinorAxisLength
    stats.Orientation
    stats.Solidity
    stats.Extent
    stats.Eccentricity
    stats.ConvexArea
    stats.FilledArea]';

Xn = [
    'center of grav i [px]   '
    'center of grav j [px]   '
    'Height [px]             '
    'Width [px]              '
    'Area [px]               '
    'Perimeter [px]          '
    'Roundness               '
    'Danielsson factor       '
    'Euler Number            '
    'Equivalent Diameter [px]'
    'MajorAxisLength [px]    '
    'MinorAxisLength [px]    '
    'Orientation [grad]      '
    'Solidity                '
    'Extent                  '
    'Eccentricity            '
    'Convex Area [px]        '
    'Filled Area [px]        '];


