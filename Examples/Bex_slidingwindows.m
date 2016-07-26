%
% Toolbox: Balu
%    Example: Welding defect detection using Sliding Windows
% 
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl



% Feature Extraction
I                = imread('testimg7.bmp'); % grayvalue image
J                = imread('testimg8.bmp'); % ideal segmentation
I1               = I(:,1:250);
J1               = J(:,1:250);
bf.name          = 'lbp';                  % definition of
bf.options.show  = 0;                      % LBP features
bf.options.vdiv  = 1;
bf.options.hdiv  = 1;
opf.b            = bf;
opf.colstr       = 'g';                    % image is grayvalue
options.opf      = opf;
options.selec    = 0;                      % all features
options.m        = 24;                     % size of a window mxm
options.n0       = 100;                    % number of 0 windows
options.n1       = 100;                    % number of 1 windows
options.th0      = 0.02;                   % threshold for 0
options.th1      = 0.02;                   % threshold for 1
options.show     = 1;
[X,d,Xn]         = Bfx_randomsliwin(I1,J1,options);

% Feature Selection
op.m             = 10;                     % 10 features will be selected
op.s             = 0.75;                   % only 75% of sample will be used
op.show          = 0;                      % no display results
op.b.name        = 'fisher';               % SFS with Fisher
selec            = Bfs_balu(X,d,op);       % selected features index
f                = X(:,selec);             % selected features

% Training
b.name           = 'lda';
b.options.p      = [];
opc              = b;
opc              = Bcl_structure(f,d,opc);

% Detection
options.opc      = opc;
options.nm       = 6;                     % shifted by 24/6=4 pixels
options.Dth      = 24;                    % 24x24 pixels
options.selec    = selec;
I2               = I(:,251:500);
[Dmap,Dbin]      = Bim_segsliwin(I2,options);
% >>> Compare with ideal detection J2 = J(:,251:500);
% >>> This is an example only, better results can be obtained using
% >>> more training images and more features.


