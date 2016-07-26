% [X,Xn] = Bfx_int(I,R,b)
%
% Toolbox: Balu
%
%    Intensity feature extraction.
%
%    This function calls intensity feature extraction procedures of image I
%    according binary image R. See example to see how it works.
%
%    X is the feature matrix (one feature per column, one sample per row),
%    Xn is the list of feature names (see Examples to see how it works).
%
%   Example 1: Extraction of one region image in grayvalue images
%      b(1).name = 'gabor';       b(1).options.show=1;         % Gabor features
%                                 b(1).options.Lgabor  = 8;    % number of rotations
%                                 b(1).options.Sgabor  = 8;    % number of dilations (scale)
%                                 b(1).options.fhgabor = 2;    % highest frequency of interest
%                                 b(1).options.flgabor = 0.1;  % lowest frequency of interest
%                                 b(1).options.Mgabor  = 21;   % mask size
%      b(2).name = 'basicint';    b(2).options.show=1;         % Basic intensity features
%      b(3).name = 'lbp';         b(3).options.show=1;         % LBP
%                                 b(3).options.vdiv = 2;       % vertical div
%                                 b(3).options.hdiv = 2;       % horizontal div
%      options.b = b;
%      options.colstr = 'i';                                   % gray image
%      I = imread('testimg1.jpg');                             % input image
%      R = Bim_segbalu(I);                                     % segmentation
%      J = rgb2gray(I);                                        % grayvalue image
%      [X,Xn] = Bfx_int(J,R,options);                          % intensity features
%      Bio_printfeatures(X,Xn)
%
%   Example 2: Extraction of multiple regions image
%      b(1).name = 'huint';       b(1).options.show=1;         % Hu moments
%      b(2).name = 'basicint';    b(2).options.show=1;         % Basic intensity features
%      b(3).name = 'lbp';         b(3).options.show=1;         % LBP
%                                 b(3).options.vdiv = 2;       % vertical div
%                                 b(3).options.hdiv = 2;       % horizontal div
%      b(4).name = 'contrast';    b(4).options.show = 1;       % Contrast
%                                 b(4).options.neighbor = 1;   % neighborhood is a window
%                                 b(4).options.param = 1.5;    % 1.5 height x 1.5 width
%      I = imread('rice.png');                                 % input image
%      [R,m] = Bim_segmowgli(I,ones(size(I)),40,1.5);          % segmentation
%      options.b = b;
%      options.colstr = 'i';                                   % gray image
%      [X,Xn] = Bfx_int(I,R,options);                          % intensity features
%      figure; hist(X(:,9));xlabel([Xn(9,:)])                  % std dev  histogramm
%      figure; hist(X(:,253));xlabel([Xn(253,:)])              % contrast Ks histogramm
%
%   Example 3: Extraction of one region image in RGB images
%      b(1).name = 'gabor';       b(1).options.show=1;         % Gabor features
%                                 b(1).options.Lgabor  = 8;    % number of rotations
%                                 b(1).options.Sgabor  = 8;    % number of dilations (scale)
%                                 b(1).options.fhgabor = 2;    % highest frequency of interest
%                                 b(1).options.flgabor = 0.1;  % lowest frequency of interest
%                                 b(1).options.Mgabor  = 21;   % mask size
%      b(2).name = 'basicint';    b(2).options.show=1;         % Basic intensity features
%      b(3).name = 'lbp';         b(3).options.show=1;         % LBP
%                                 b(3).options.vdiv = 2;       % vertical div
%                                 b(3).options.hdiv = 2;       % horizontal div
%      options.b = b;
%      options.colstr = 'rgb';                                 % R image
%      I = imread('testimg1.jpg');                             % input image
%      R = Bim_segbalu(I);                                     % segmentation
%      [X,Xn] = Bfx_int(I,R,options);                          % intensity features
%      Bio_printfeatures(X,Xn)
%
%   See also Bfx_basicint, Bfx_haralick, Bfx_gabor, Bfx_dct, Bfx_fourier,
%            Bfx_huint, Bfx_lbp, Bfx_contrast, Bfx_clp, Bfx_phog,
%            Bfx_files.
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [X,Xn] = Bfx_int(I,R,options)




I = double(I);

[N,M,P] = size(I);

c = size(I,3);


if nargin==2;
    options = R;
    R = ones(size(I));
end

if isempty(R)
    R = ones(N,M);
end

b    = options.b;
colstr = options.colstr;

n  = length(b);
m  = int16(max(R(:)));
X  = [];
Xn = [];
for j=1:m
    Rj = R==j;
    Xj = [];
    if j==1
        Xn = [];
    end
    for i=1:n
        s = b(i).name;
        if sum(s(1:2)=='Bf')~=2
            s = ['Bfx_' s];
        end
        if  ~exist(s,'file')
            error(sprintf('Bfx_extraction: function %s does not exist.',b(i).name))
        end
        % if or(compare(s,'Bfx_lbphogi')==0,compare(s,'Bfx_hogi')==0)
        % if (compare(s(1:8),'Bfx_lbpi')*compare(s,'Bfx_lbphogi')*compare(s,'Bfx_hogi'))==0
        if (compare(s,'Bfx_lbpi')*compare(s,'Bfx_lbphogi')*compare(s,'Bfx_hogi'))==0
            c = 1;
            ifull = 1;
        else
            c = P;
            ifull = 0;
        end
        for k=1:c
            %        for i=1:n
            if ifull
                % tic
                [Xk,Xnk] = feval(s,I,Rj,b(i).options);
                % toc
            else
                [Xk,Xnk] = feval(s,I(:,:,k),Rj,b(i).options);
            end
            Xj = [Xj Xk];
            if j==1
                nk = length(Xk);
                Xns = [ones(nk,1)*[colstr(k) '-'] Xnk ];
                Xn = [Xn; Xns(:,1:24)];
            end
        end
    end
    X = [X;Xj];
end