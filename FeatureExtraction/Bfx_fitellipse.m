% [X,Xn] = Bfx_fitellipse(R,options)
% [X,Xn] = Bfx_fitellipse(R)
%
% Toolbox: Balu
%
%    Fit ellipse for the boundary of a binary image R.
%
%    options.show = 1 display mesagges.
%
%    X is a 6 elements vector:
%      X(1): Ellipse-centre i direction
%      X(2): Ellipse-centre j direction
%      X(3): Ellipse-minor axis
%      X(4): Ellipse-major axis
%      X(5): Ellipse-orientation
%      X(6): Ellipse-eccentricity
%      X(7): Ellipse-area
%
%      Xn is the list of feature names.
%
%      Xn is the list of feature names.
%
%   Reference:
%      Fitzgibbon, A.; Pilu, M. & Fisher, R.B. (1999): Direct Least Square
%      Fitting Ellipses, IEEE Trans. Pattern Analysis and Machine
%      Intelligence, 21(5): 476-480.
%
%   Example:
%      I = imread('testimg1.jpg');            % input image
%      R = Bim_segbalu(I);                    % segmentation
%      [X,Xn] = Bfx_fitellipse(R);            % ellipse features
%      Bio_printfeatures(X,Xn)
%
%   See also Bfx_basicgeo, Bfx_hugeo, Bfx_gupta, Bfx_flusser.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function  [X,Xn] = Bfx_fitellipse(R,options)

E = bwperim(R,4);
[Y,X] = find(E==1);               % pixel of perimeter in (i,j)
if length(X)>5

    if ~exist('options','var')
        options.show = 0;
    end

    if options.show == 1
        disp('--- extracting ellipse features...');
    end

    % normalize data
    mx = mean(X);
    my = mean(Y);
    sx = (max(X)-min(X))/2;
    sy = (max(Y)-min(Y))/2;
    x = (X-mx)/sx;
    y = (Y-my)/sy;

    % Build design matrix
    D = [ x.*x  x.*y  y.*y  x  y  ones(size(x)) ];

    [U,S,V] = svd(D);
    A = V(:,6);

    % unnormalize
    a = [
        A(1)*sy*sy,   ...
        A(2)*sx*sy,   ...
        A(3)*sx*sx,   ...
        -2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy,   ...
        -A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy,   ...
        A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my   ...
        - A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my   ...
        + A(6)*sx*sx*sy*sy   ...
        ]';

    a = a/a(6);

    % get ellipse orientation
    alpha = atan2(a(2),a(1)-a(3))/2;

    % get scaled major/minor axes
    ct = cos(alpha);
    st = sin(alpha);
    ap = a(1)*ct*ct + a(2)*ct*st + a(3)*st*st;
    cp = a(1)*st*st - a(2)*ct*st + a(3)*ct*ct;

    % get translations
    T = [[a(1) a(2)/2]' [a(2)/2 a(3)]'];
    mc = -inv(2*T)*[a(4) a(5)]';

    % get scale factor
    val = mc'*T*mc;
    scale = abs(1 / (val- a(6)));

    % get major/minor axis radii
    ae  = 1/sqrt(scale*abs(ap));
    be  = 1/sqrt(scale*abs(cp));
    ecc = ae/be; % eccentricity
    ar  = pi*ae*be;

    X = [ mc' ae be alpha ecc ar];
else
    X = [0 0 0 0 0 0 0];
    disp('Warning: Bfx_fitellipse does not have enough points to fit');
end
Xn = [
    'Ellipse-centre i [px]   '
    'Ellipse-centre j [px]   '
    'Ellipse-minor ax [px]   '
    'Ellipse-major ax [px]   '
    'Ellipse-orient [rad]    '
    'Ellipse-eccentricity    '
    'Ellipse-area [px]       '
    ];


