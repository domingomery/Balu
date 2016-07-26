% [X,d,Xn,x] = Bfx_randomsliwin(I,J,options)
%
% Toolbox: Balu
%
%    Feature extraction of random sliding windows.
%    This program select automatically detection windows sized mxm
%    with label '1' and lable '0'. For each window
%    Balu intensity features are extracted.
%
%    Input:
%    I             original image (more than one channel is allowed)
%    J             ideal segmentation
%    options.opf   feature extraction options (see example)
%    options.selec selected features, selec = 0 means all features
%    options.m     sliding window size in pixels (mxm)
%    options.n0    number of '0' windows
%    options.n1    number of '1' windows
%    options.ROI   region of interest where the windows are extracted
%    options.th0   if the number of '1' in detection window/m^2 < th0 a '0'
%                  sample is selected
%    options.th1   if the number of '1' in detection window/m^2 >=th1 a '1'
%                  sample is selected
%    options.show  display detected windows
%
%    Output:
%    X     feature values
%    Xn    feature names
%    d     ideal classification (0 or 1) of each sample
%    x     ceter of mass (j,i) of each patch
%
%    Example:
%      I1                 = imread('testimg7.bmp');  % grayvalue image
%      [I_on,I2]          = Bim_cssalient(I1,1,0);   % saliency map
%      I(:,:,1)           = I1;                      % channel 1
%      I(:,:,2)           = I2;                      % cahnnel 2
%      J                  = imread('testimg8.bmp');  % ideal segmentation
%      bf(1).name         = 'lbp';                   % definition of
%      bf(1).options.show = 0;                       % first features
%      bf(1).options.vdiv = 1;
%      bf(1).options.hdiv = 1;
%      bf(2).name         = 'basicint';              % definition of
%      bf(2).options.show = 0;                       % second features
%      bf(2).options.mask = 5;
%      opf.b              = bf;
%      opf.colstr         = 'gs';                    % chn 1,2 are gray,sal
%      options.opf        = opf;
%      options.selec      = 0;                       % all features
%      options.m          = 24;                      % size of a window mxm
%      options.n0         = 100;                     % number of 0 windows
%      options.n1         = 100;                     % number of 1 windows
%      options.th0        = 0.02;                    % threshold for 0
%      options.th1        = 0.02;                    % threshold for 1
%      options.show       = 1;
%      [X,d,Xn]           = Bfx_randomsliwin(I,J,options);
%
%  See also Bim_segsliwin.
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [X,d,Xn,x] = Bfx_randomsliwin(I,J,options)

warning off

opf   = options.opf;
% selec = options.selec;
m     = options.m;     % size of a window mxm
n0    = options.n0;
n1    = options.n1;
th0   = options.th0;
th1   = options.th1;
show  = options.show;
if isfield(options,'win');
    win   = options.win;
else
    win   = 30;
end

if isfield(options,'roi');
    ROI   = options.roi;
else
    ROI = ones(size(I));
end

if isfield(options,'selec');
    selec  = options.selec;
else
    selec = 0;
end



[N,M]=size(I(:,:,1));

if show==1
    close all
    figure(1)
    imshow(I(:,:,1),[]);title('Original');
    pause(10)
    hold on
    %figure(2)
    %imshow(J);title('Real Defects');
    %figure(1)
end


m1  = m-1;
md  = (m1-1)/2;
m2  = m*m;
R   = ones(m,m);

if show~=-1
    ff = Bio_statusbar('Extracting patches');
end

ft   = Bfx_int(I(1:win,1:win),opf);
nf   = size(ft,2);

% f   = [];
% x   = [];
nn = n0+n1;
f = zeros(nn,nf);
x = zeros(nn,2);
d = [zeros(n0,1); ones(n1,1)];

k = 0;


if n0>0
    % Extracting features for '0' detection windows
    if show==1
        disp('Extracting features for 0 detection windows...');
    end
    i   = 0;
    while i<n0
        i1 = fix(rand*(N-m1))+1;
        j1 = fix(rand*(M-m1))+1;
        rj = sum2(ROI(i1:i1+m1,j1:j1+m1))/m2;
        if rj>0.90
            wj = J(i1:i1+m1,j1:j1+m1);
            t = sum(wj(:))/m2;
            if t<th0
                wij  = I(i1:i1+m1,j1:j1+m1,:);
                [ft,fn]   = Bfx_int(wij,R,opf);
                % f = [f;ft 0];
                % x = [x;j1+md i1+md];
                k = k+1;
                f(k,:) = ft;
                x(k,:) = [j1+md i1+md];
                i = i+1;
                if show~=-1
                    ff = Bio_statusbar(k/nn,ff);
                end
                if show==1
                    % fprintf('0: %d/%d\n',i,n0);
                    plot([j1 j1 j1+m1 j1+m1 j1],[i1 i1+m1 i1+m1 i1 i1],'g')
                    drawnow
                end
            end
        end
    end
end

if n1>0
    % Extracting features for '1' detection windows
    if show==1
        disp('Extracting features for 1 detection windows...');
    end
    i = 0;
    while i<n1
        i1 = fix(rand*(N-m1))+1;
        j1 = fix(rand*(M-m1))+1;
        rj = sum2(ROI(i1:i1+m1,j1:j1+m1))/m2;
        if rj>0.95
            wj = J(i1:i1+m1,j1:j1+m1);
            t = sum(wj(:))/m2;
            if t>=th1
                wij     = I(i1:i1+m1,j1:j1+m1,:);
                [ft,fn] = Bfx_int(wij,R,opf);
                % f       = [f;ft 1];
                % x       = [x;j1+md i1+md];
                k = k+1;
                f(k,:) = ft;
                x(k,:) = [j1+md i1+md];
                i       = i+1;
                if show~=-1
                    ff = Bio_statusbar(k/nn,ff);
                end
                
                if show==1
                    % fprintf('1: %d/%d\n',i,n1);
                    plot([j1 j1 j1+m1 j1+m1 j1],[i1 i1+m1 i1+m1 i1 i1],'r')
                    drawnow
                end
            end
        end
    end
end
if sum(selec)==0
    X  = f;
    Xn = fn;
else
    X  = f(:,selec);
    Xn = fn(selec,:);
end
if show~=-1
    delete(ff);
end



