% DELETE ALL
clt

% IMAGE AQUISITION DEFINITION
f.path          = '/Users/domingomery/matlab/espinas2';  % current directory or a path directory
f.prefix        =  '*';   
f.extension   =  'png';
f.imgmin        = 1;
f.imgmax        = 6;


% GEOMETRIC FEATURE EXTRACTION DEFINTION
k = 0;

k = k+1;
b(k).name = 'basicgeo';                     % basic geometric features
             b(k).options.show        = 1;    % display results

k = k+1;
b(k).name = 'fitellipse';                     % elliptic features
             b(k).options.show        = 1;    % display results

k = k+1;
b(k).name = 'fourierdes';                     % Fourier descriptors
             b(k).options.Nfourierdes = 16;   % number of descriptors
             b(k).options.show        = 1;    % display results

k = k+1;
b(k).name = 'hugeo';                          % Hu moments
             b(k).options.show        = 1;    % display results

k = k+1;
b(k).name = 'flusser';                        % Flusser moments
             b(k).options.show        = 1;    % display results

k = k+1;
b(k).name = 'gupta';                          % Gupta moments
             b(k).options.show        = 1;    % display results


opg.b            = b;
opi.segmentation = 'Bim_segbalu';
opi.param        = 0;
             
             
             
             
% INTENSITY FEATURE EXTRACTION DEFINTION
k = 0;

k = k+1;
b(k).name = 'basicint';                       % basic intensity features
             b(k).options.show        = 1;    % display results

k = k+1;
b(k).name = 'contrast';                     % contrast features
             b(k).options.neighbor    = 2;    % neigborhood is imdilate
             b(k).options.param       = 5;    % with 5x5 mask
             b(k).options.show        = 1;    % display results

k = k+1;
b(k).name = 'haralick';                       % statistical texture features
             b(k).options.dharalick   = 1:5;  % for 1, 2, ... 5 pixels
             b(k).options.show        = 1;    % display results

k = k+1;
b(k).name = 'lbp';                            % local binary batterns (LBP)
             b(k).options.vdiv        = 1;    % one vertical divition
             b(k).options.hdiv        = 1;    % one horizontal divition
             b(k).options.samples     = 8;    % number of neighbor samples
             b(k).options.show        = 1;    % display results

k = k+1;
b(k).name = 'lbp';                            % semantic LBP
             b(k).options.vdiv        = 1;    % one vertical divition
             b(k).options.hdiv        = 1;    % one horizontal divition
             b(k).options.semantic    = 1;    % semantic LBP
             b(k).options.samples     = 8;    % number of neighbor samples
             b(k).options.sk          = 0.5;  % angle sampling
             b(k).options.show        = 1;    % display results

k = k+1;
b(k).name = 'dct';                            % Discrete Cosinus Transform
             b(k).options.Ndct        = 64;   % imresize vertical
             b(k).options.Mdct        = 64;   % imresize horizontal
             b(k).options.mdct        = 4;    % imresize frequency vertical
             b(k).options.ndct        = 4;    % imresize frequency horizontal
             b(k).options.show        = 1;    % display results
            
k = k+1;
b(k).name = 'fourier';                        % Discrete Fourier Transform
             b(k).options.Nfourier    = 64;   % imresize vertical
             b(k).options.Mfourier    = 64;   % imresize horizontal
             b(k).options.mfourier    = 4;    % imresize frequency vertical
             b(k).options.nfourier    = 4;    % imresize frequency horizontal
             b(k).options.show        = 1;    % display results
       
       
k = k+1;
b(k).name = 'gabor';                          % Gabor features
             b(k).options.Lgabor      = 8;    % number of rotations
             b(k).options.Sgabor      = 8;    % number of dilations (scale)
             b(k).options.fhgabor     = 2;    % highest frequency of interest
             b(k).options.flgabor     = 0.1;  % lowest frequency of interest
             b(k).options.Mgabor      = 21;   % mask size
             b(k).options.show        = 1;    % display results

k = k+1;
b(k).name = 'huint';                          % Hu-moments with intensity  
             b(k).options.show        = 1;    % display results

opi.b = b;
opi.colstr = 'g';

             

[X,Xn,S] = Bfx_files(f,opg,opi);