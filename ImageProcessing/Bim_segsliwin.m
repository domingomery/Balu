% [Dmap,Dbin] = Bim_segsliwin(I,options)
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
%    options.opf   feature extraction options (see example)
%    options.opc   classifier options (see example)
%    options.selec selected features, selec = 0 means all features
%    options.m     sliding window size in pixels (mxm)
%    options.nm    shifted portion (shifted pixels = m/nm)
%    options.Dth   Threshold (if Dth==0 then Dth = (nm^2)/2-1)
%    options.show  display detected windows
%
%    Output:
%    Dbin          Binary image for the detection
%    Dmap          Detection map (Dbin = Dmap>options.Dth).
%
%    Example:
%      % Feature Extraction
%      I                = imread('testimg7.bmp'); % grayvalue image
%      J                = imread('testimg8.bmp'); % ideal segmentation
%      I1               = I(:,1:250);
%      J1               = J(:,1:250);
%      bf.name          = 'lbp';                  % definition of
%      bf.options.show  = 0;                      % LBP features
%      bf.options.vdiv  = 1;
%      bf.options.hdiv  = 1;
%      bf.name          = 'basicint';             % definition of
%      bf.options.show  = 0;                      % Basic intensity features
%      bf.options.mask  = 5;
%      opf.b            = bf;
%      opf.colstr       = 'g';                    % image is grayvalue
%      options.opf      = opf;
%      options.selec    = 0;                      % all features
%      options.m        = 24;                     % size of a window mxm
%      options.n0       = 100;                    % number of 0 windows
%      options.n1       = 100;                    % number of 1 windows
%      options.th0      = 0.02;                   % threshold for 0
%      options.th1      = 0.02;                   % threshold for 1
%      options.show     = 1;
%      [X,d,Xn]         = Bfx_randomsliwin(I1,J1,options);
%
%      % Feature Selection
%      op.m             = 3;                     % 10 features will be selected
%      op.s             = 0.75;                   % only 75% of sample will be used
%      op.show          = 0;                      % no display results
%      op.b.name        = 'fisher';               % SFS with Fisher
%      selec            = Bfs_balu(X,d,op);       % selected features index
%      fx               = X(:,selec);             % selected features
%
%      % Training
%      b.name           = 'lda';
%      b.options.p      = [];
%      opc              = b;
%      opc              = Bcl_structure(fx,d,opc);
%
%      % Detection
%      options.opc      = opc;
%      options.nm       = 6;                     % shifted by 24/6=4 pixels
%      options.Dth      = 24;                    % 24x24 pixels
%      options.selec    = selec;
%      I2               = I(:,251:500);
%      [Dmap,Dbin]      = Bim_segsliwin(I2,options);
%      % >>> Compare with ideal detection J2 = J(:,251:500);
%      % >>> This is an example only, better results can be obtained using
%      % >>> more training images and more features.
%
%
% (c) D.Mery, PUC-DCC, 2011-2014
% http://dmery.ing.puc.cl

function [Dmap,Dbin,MMM] = Bim_segsliwin(I,options)

opf   = options.opf;
opc   = options.opc;
selec = options.selec;
m     = options.m;
nm    = options.nm;
show  = options.show;
Dth   = options.Dth;

if isfield(options,'roi');
    ROI = options.roi;
else
    ROI = ones(size(I));
end

if isfield(options,'detection');
    dtec = options.detection;
else
    dtec = 1;
end


N = size(I,1);
M = size(I,2);

if show==1
    close all
    figure(1)
    imshow(I(:,:,1),[]);
    title('Original');
    drawnow
end


% Parameters
if Dth==0
    Dth   = nm*nm/2-1;
end
dm   = round(m/nm);            % overlapping
m1   = m-1;
m2   = m*m;
R    = ones(m,m);


% Feature Extraction of Sliding Windows
si1 = 1:dm:N-m;
sj1 = 1:dm:M-m;

wij = I(1:m,1:m,:);
ft  = Bfx_int(wij,R,opf);
nft = length(si1)*length(sj1);
mft = length(ft);
fj  = zeros(nft,mft);
k1  = 1;
ff = Bio_statusbar('Sliding Windows');
for i1=1:dm:N-m
    ff = Bio_statusbar(i1/N,ff);
    if show
        fprintf('processing image row %d/%d...\n',i1,N-m)
    end
    for j1 = 1:dm:M-m
        rj = sum2(ROI(i1:i1+m1,j1:j1+m1))/m2;
        if rj>0.95
            fj(k1,:) = Bfx_int(I(i1:i1+m1,j1:j1+m1,:),R,opf);
            k1       = k1+1;
        end
    end
end
delete(ff);

if sum(selec)==0
    X = fj;
else
    X = fj(:,selec);
end

% Classification of all Sliding Windows
dt      = Bcl_structure(X,opc);

% Detection Map
j = 0;
Dmap = zeros(N,M);
ff = Bio_statusbar('Detecting');
for i1=1:dm:N-m
    ff = Bio_statusbar(i1/N,ff);
    for j1 = 1:dm:M-m
        rj = sum2(ROI(i1:i1+m1,j1:j1+m1))/m2;
        if rj>0.95
            j = j+1;
            if dt(j)==dtec
                Dmap(i1:i1+m1,j1:j1+m1)=Dmap(i1:i1+m1,j1:j1+m1)+1;
            end
        end
    end
end
delete(ff);

% Binary Detection
Dbin = Dmap>Dth;

% Display Results
tini = 0;
jjj = 1;
if show
    figure(2)
    imshow(I(:,:,1),[]);
    title('Detected Sliding Windows');
    hold on
    j=0;
    for i1=1:dm:N-m
        for j1 = 1:dm:M-m
            rj = sum2(ROI(i1:i1+m1,j1:j1+m1))/m2;
            if rj>0.95
                j = j+1;
                if dt(j)==dtec
                    x = [j1 j1+m1 j1+m1 j1 j1];
                    y = [i1 i1 i1+m1 i1+m1 i1];
                    plot(x,y,'r')
                end
            end
        end
        if tini==0
            %for kkkk=1:3
            %    kkkk
            %    pause(1)
            %end
            tini=1;
        end
        pause(0.01)
        MMM(jjj) = getframe;
        jjj=jjj+1;
    end
    % enterpause

    figure(3)
    % imshow(Dmap,[]);
    Dmap = log(abs(Dmap)+1);
    imshow(uint8(double(Dmap)/max2(Dmap)*63),jet);
    title('Detection Map');

    figure(4)
    drawnow
    hold on
    E = bwperim(Dbin);
    Bio_edgeview(I(:,:,1),imdilate(E,ones(3,3)))
    title('Binary Detection');
end
