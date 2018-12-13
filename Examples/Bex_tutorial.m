% Tutorial for Balu
% Oficial webpage: http://dmery.ing.puc.cl/index.php/balu
%
% (c) Domingo Mery
% Department of Computer Science
% Universidad Catolica de Chile
% http://dmery.ing.puc.cl
%
% September 2015
%
%
% 0: Download
% >>>> http://dmery.ing.puc.cl/index.php/balu/download/
% See slides 2015-BaluTutorial.pptx for modules
% https://www.dropbox.com/s/83taze5ecq5wsfb/2015-Balu-Tutorial.pptx?dl=0
%

%% 1.1: IMAGE PROCESSING - SALIENCY
%
%  [J_on, J_off] = Bim_cssalient( I, type, show )
%
%     compute center surround saliency image from and gray scale image follow
%     mehtod mentioned in Montabone & Soto (2010). Method precompute an
%     integral image to do all calculations.
%
%     Inputs:
%     I:  gray scale image. If depth of I is higher that 1 the algoritm
%         transform the image into gray scale version.
%     type:   scales method
%     show:   plot options. 1 shows both saliency images. 0 shows anything.
%
%     Output:
%     J_off: Surround-Center
%     J_on: Center-Surround
%
clt
I = imread('testimg7.bmp');
[J_on,J_off] = Bim_cssalient(I,1,1);

%% 1.2: IMAGE PROCESSING - SEGMENTATION
%
% [R,E,J] = Bim_segbalu(I,p)
%    Segmentation of an object with homogeneous background.
%
%    I: input image
%    p: threshold (default: p=-0.05) with p between -1 and 1.
%        A positive value is used to dilate the segmentation,
%        the negative to erode.
%    R: binary image of the object
%    E: binary image of the edge of the object
%    J: high contrast image of I.
%
clt
I = imread('testimg1.jpg');
R = Bim_segbalu(I);
figure(1)
imshow(I); title('test image')
figure(2)
imshow(R); title('segmented image')
enterpause
% Segmentation and display results
I = imread('testimg2.jpg');
Bio_segshow(I,'Bim_segbalu')

%% 1.3: IMAGE PROCESSING - EQUALIZATION
%
%   Y = Bim_equalization(X,show)
%
%      Enhancement of a grayvalue image forcing an uniform histogram.
%
%      Input data:
%         X grayvalue image.
%
%      Output:
%         Y: enhanced image so that histogram of Y is perfectlty uniformed
%         distributed. Y is uint8.
%
clt
X = imread('tire.tif');
Y = Bim_equalization(X,1);
figure(3)
imshow([X Y])

%% 1.4: IMAGE PROCESSING - SPOT SEGMENTATION
% [F,m] = Bim_segmowgli(J,R,Amin,sig)
%
%    Segmentation of regions in image J using LoG edge detection.
%    R   : binary image of same size of J that indicates the piexels where
%          the segmentation will be performed. Default R = ones(size(J));
%    Amin: minimum area of the segmented details.
%    sig : sigma of LoG edge detector.
%    F   : labeled image of the segmentation.
%    m   : numbers of segmented regions.
%
%    Example 1:
clt
I = imread('rice.png');
figure(1);imshow(I);title('test image');
F = Bim_segmowgli(I,[],40,1.5);
figure(2);Bio_edgeview(I,F>0);title('segmented image');
enterpause

%   Example 2:
I = imread('testimg4.jpg');
figure(1);imshow(I);title('test image');
R = Bim_segbalu(I,-0.1);
figure(2);imshow(R);title('segmented object');
G = I(:,:,2);
[F,m] = Bim_segmowgli(G,R,30,2);
figure(3);imshow(F,[]);title('segmented regions')

%% 1.5: IMAGE PROCESSING - LOOK UP TABLE
%  function Y = LUT(X,T,show)
%
%      Look Up Table transformation for grayvalue images.
%
%      Input data:
%         X grayvalue image.
%         T look uo table
%         show display results
%
%      Output:
%         Y transformed image
%
clt
load clown
T = 256*ones(256,1);
T(1:81) = (1:81).*(1:81)/81^2*255; % look up table
Y = Bim_LUT(X,T,1);

%% 1.6: IMAGE PROCESSING - IMAGE RESTORATION
%
% Image restoration using MINIO criterium.
%
%    G : blurred image
%    h : PSF
%    method = 1: minio ||f_N-g|| -> min (default)
%           = 2: ||f|| -> min
%    Fs: restored image
clt
F = double(imread('circuit.tif'));
n = 45;
h = ones(1,n)/n;
G = conv2(F,h,'valid');
Fs = Bim_resminio(G,h);
figure(1)
imshow(F,[]); title('original image')
figure(2)
imshow(G,[]); title('blurred image')
figure(3)
imshow(Fs,[]); title('restored image')

%% 1.7: IMAGE PROCESSING - COLOR SEGMENTATION USING K-MEANS
%
% [R,J] = Bim_segkmeans(I,k,r,show)
%
%      Segmentation of color images using kmeans.
%
%   I: input image
%   k: number of clusters
%   r: resize image
%   show: 1 means intermediate results will be displayed
%
clt
I = imread('testimg9.jpg');
[R,J] = Bim_segkmeans(I,3,1,1);


%% 2.1: FEATURE EXTRACTION - BASIC GEOMETRIC FEATURES
%
%[X,Xn] = Bfx_geobasic(R,options)
%
% Standard geometric features of a binary image R. This function calls
%      regionprops of Image Processing Toolbox.
%
%      options.show = 1 display mesagges.
%
%      X is the feature vector
%      Xn is the list of feature names.
clt
I = imread('testimg1.jpg');            % input image
[R,E] = Bim_segbalu(I);                % segmentation
Bio_edgeview(I,E,[0 1 0],3)
[X,Xn] = Bfx_basicgeo(R);              % basic geometric features
Bio_printfeatures(X,Xn)

%% 2.2: FEATURE EXTRACTION - FOURIER DESCRIPTORS
%   function [X,Xn] = Bfx_fourierdes(R,options)
%
%      Computes the Fourier descriptors of a binary image R.
%
%      options.show = 1 display mesagges.
%      options.Nfourierdes number of descriptors.
%
%      X is the feature vector
%      Xn is the list of feature names.
%
I = imread('testimg1.jpg');            % input image
[R,E] = Bim_segbalu(I);                % segmentation
Bio_edgeview(I,E,[0 1 0],3)
op.show = 1;
op.Nfourierdes = 16;                   % Number of Fourier Descriptors
[X,Xn] = Bfx_fourierdes(R,op);         % Fourier descriptors
Bio_printfeatures(X,Xn)

%% 2.3: FEATURE EXTRACTION - CROSSING LINE PROFILE
%  [X,Xn] = Bfx_clp(I,R,options)
%  [X,Xn] = Bfx_clp(I,options)
%
%      Crossing Line Profile.
%
%      X is the features vector, Xn is the list of feature names(see Example
%      to see how it works).
clt
options.show    = 1;                    % display results
options.ng      = 32;                   % windows resize
I = imread('testimg4.jpg');             % input image
J = I(395:425,415:442,1);               % region of interest (red)
R = J>135;                              % segmentation
figure;imshow(J,[])
figure;imshow(R)
[X,Xn] = Bfx_clp(J,R,options);          % CLP features
Bio_printfeatures(X,Xn)

%% 2.4: FEATURE EXTRACTION - GABOR
%   [X,Xn] = Bfx_gabor(I,R,options)
%   [X,Xn] = Bfx_gabor(I,options)
%
%      Gabor features
%
%      X is the features vector, Xn is the list of feature names (see Example
%      to see how it works).
clt
options.Lgabor  = 8;                 % number of rotations
options.Sgabor  = 8;                 % number of dilations (scale)
options.fhgabor = 2;                 % highest frequency of interest
options.flgabor = 0.1;               % lowest frequency of interest
options.Mgabor  = 21;                % mask size
options.show    = 1;                 % display results
I = imread('testimg1.jpg');          % input image
[R,E] = Bim_segbalu(I);              % segmentation
Bio_edgeview(I,E,[0 1 0],3)
J = I(:,:,2);                        % green channel
[X,Xn] = Bfx_gabor(J,R,options);     % Gabor features
Bio_printfeatures(X,Xn)

%% 2.5: FEATURE EXTRACTION - LBP
%   [X,Xn] = Bfx_lbp(I,R,options)
%   [X,Xn] = Bfx_lbp(I,options)
%
%      Local Binary Patterns features
%
%      X is the features vector, Xn is the list of feature names (see Example
%      to see how it works).
%
%      It calculates the LBP over the a regular grid of patches. The function
%      uses Heikkila & Ahonen (see http://www.cse.oulu.fi/MVG/Research/LBP).
%
%      It returns a matrix of uniform lbp82 descriptors for I, made by
%      concatenating histograms of each grid cell in the image.
%      Grid size is options.hdiv * options.vdiv
%
%      R is a binary image or empty. If R is given the lbp will be computed
%      the corresponding pixles R==0 in image I will be set to 0.
%
%       Output:
%       X is a matrix of size ((hdiv*vdiv) x 59), each row has a
%           histogram corresponding to a grid cell. We use 59 bins.
%       options.x of size hdiv*vdiv is the x coordinates of center of ith grid cell
%       options.y of size hdiv*vdiv is the y coordinates of center of ith grid cell
%       Both coordinates are calculated as if image was a square of side length 1.
clt
options.vdiv = 1;                  % one vertical divition
options.hdiv = 1;                  % one horizontal divition
options.semantic = 0;              % classic LBP
options.samples  = 8;              % number of neighbor samples
options.mappingtype = 'u2';        % uniform LBP
I = imread('testimg1.jpg');        % input image
J = I(120:219,120:239,2);          % region of interest (green)
figure(1);imshow(J,[])             % image to be analyzed
[X,Xn] = Bfx_lbp(J,[],options);    % LBP features
figure(2);bar(X)                   % histogram

%% 2.6 FEATURE EXTRACTION - SEVERAL GEOMETRIC FEATURES
%
% In this example we show how to extract several features of many regions.
clt

% Input Image
I = imread('rice.png');
% Segmentation
[R,m] = Bim_segmowgli(I,ones(size(I)),40,1.5);


% Definition of features to be extracted
b(1).name = 'hugeo';       b(1).options.show=1;                                 % Hu moments
b(2).name = 'basicgeo';    b(2).options.show=1;                                 % basic geometric fetaures
b(3).name = 'fourierdes';  b(3).options.show=1; b(3).options.Nfourierdes = 16;  % Fourier Descriptors
options.b = b;

% Feature extraction
[X,Xn] = Bfx_geo(R,options);

% Processing after feature extraction
figure; hist(X(:,12));xlabel(Xn(12,:))            % area histogram
ii = find(abs(X(:,20))<15);                       % rice orientation
K = zeros(size(R));                               % between -15 and 15 grad
for i=1:length(ii);K=or(K,R==ii(i));end
E = bwperim(K);
figure; Bio_edgeview(I,E,[1 0 1]);title('abs(orientation)<15 grad')
h = X(:,41-14:41);
figure
mesh(h)
title('Fourier Descriptors');


%% 2.7 FEATURE EXTRACTION - SEVERAL INTENSITY FEATURES
%
% [X,Xn] = Bfx_int(I,R,b)
%
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
clt
b(1).name = 'gabor';       b(1).options.show=1;         % Gabor features
b(1).options.Lgabor  = 8;    % number of rotations
b(1).options.Sgabor  = 8;    % number of dilations (scale)
b(1).options.fhgabor = 2;    % highest frequency of interest
b(1).options.flgabor = 0.1;  % lowest frequency of interest
b(1).options.Mgabor  = 21;   % mask size
b(2).name = 'basicint';    b(2).options.show=1;         % Basic intensity features
b(3).name = 'lbp';         b(3).options.show=1;         % LBP
b(3).options.vdiv = 2;       % vertical div
b(3).options.hdiv = 2;       % horizontal div
options.b = b;
options.colstr = 'i';                                   % gray image
I = imread('testimg1.jpg');                             % input image
R = Bim_segbalu(I);                                     % segmentation
J = rgb2gray(I);                                        % grayvalue image
[X,Xn] = Bfx_int(J,R,options);                          % intensity features
Bio_printfeatures(X,Xn)
enterpause

%   Example 2: Extraction of multiple regions image
clt
b(1).name = 'huint';       b(1).options.show=1;         % Hu moments
b(2).name = 'basicint';    b(2).options.show=1;         % Basic intensity features
b(3).name = 'lbp';         b(3).options.show=1;         % LBP
b(3).options.vdiv = 2;       % vertical div
b(3).options.hdiv = 2;       % horizontal div
b(4).name = 'contrast';    b(4).options.show = 1;       % Contrast
b(4).options.neighbor = 1;   % neighborhood is a window
b(4).options.param = 1.5;    % 1.5 height x 1.5 width
I = imread('rice.png');                                 % input image
[R,m] = Bim_segmowgli(I,ones(size(I)),40,1.5);          % segmentation
options.b = b;
options.colstr = 'i';                                   % gray image
[X,Xn] = Bfx_int(I,R,options);                          % intensity features
figure; hist(X(:,9));xlabel([Xn(9,:)])                  % std dev  histogramm
figure; hist(X(:,253));xlabel([Xn(253,:)])              % contrast Ks histogramm
enterpause

%   Example 3: Extraction of one region image in RGB images
clt
b(1).name = 'gabor';       b(1).options.show=1;         % Gabor features
b(1).options.Lgabor  = 8;    % number of rotations
b(1).options.Sgabor  = 8;    % number of dilations (scale)
b(1).options.fhgabor = 2;    % highest frequency of interest
b(1).options.flgabor = 0.1;  % lowest frequency of interest
b(1).options.Mgabor  = 21;   % mask size
b(2).name = 'basicint';    b(2).options.show=1;         % Basic intensity features
b(3).name = 'lbp';         b(3).options.show=1;         % LBP
b(3).options.vdiv = 2;       % vertical div
b(3).options.hdiv = 2;       % horizontal div
options.b = b;
options.colstr = 'rgb';                                 % R image
I = imread('testimg1.jpg');                             % input image
R = Bim_segbalu(I);                                     % segmentation
[X,Xn] = Bfx_int(I,R,options);                          % intensity features
Bio_printfeatures(X,Xn)




%% 2.8 FEATURE EXTRACTION USING GUI
%
clt
d = pwd;
st = [pwd '/orl/'];
orl_url = 'http://www.cl.cam.ac.uk/Research/DTG/attarchive/pub/data/att_faces.zip' ;
disp('Balu demo requires a database to test. In this case we will use');
disp('ORL database from:');
disp(' ');
disp('F.S. Samaria and A.C. Harter: Parameterisation of a stochastic model for human face');
disp('identification. In Proceedings of the Second IEEE Workshop on Applications of');
disp('Computer Vision, pages 138-142, 1994.');
disp(' ');
disp('By downloading you agree to acknowledge the source of the images by including ');
disp('the citation in your publications.');
disp(' ')
yn = input('Do you want to download the image database [ 1=yes, 0=no ] (first time answer yes)? ');
if yn==1
    fprintf('Downloading ORL data to ''%s''. This will take a while...', st);
    unzip(orl_url,'orl');
    disp(' ');
    disp('ORL download succefully!');
    disp(' ');
end

yn = input('Do you want to copy the image faces into the database [ 1=yes, 0=no ] (first time answer yes)? ');
if yn==1
    
    disp('Copying faces and croping images...');
    ft = Bio_statusbar('copying files');
    n = 40;
    for y=1:n
        ft = Bio_statusbar(y/n,ft);
        for t=1:10
            s1 = ['orl/s' num2str(y) '/' num2str(t) '.pgm' ];
            J = imresize(imread(s1),[110 90]);
            sf = ['orl/face_' num2fixstr(y,3) '_' num2fixstr(t,5) '.png' ];
            imwrite(J,sf,'png');
        end
    end
    delete(ft)
end


cd orl
disp(' ');
disp(' ');
disp('In the GUI please select png files, LBP features, Matlab output and Gray Components... and GO!');
enterpause

Bfx_gui
cd ..

%% 3.1 FEATURE SELECTION - SEQUENTIAL FORWARD SELECTION
% selec = Bfs_sfs(X,d,options)
%
%      Sequential Forward Selection for fatures X according to ideal
%      classification d. optins.m features will be selected.
%      options.b.method = 'fisher' uses Fisher objetctive function.
%      options.b.method = 'sp100' uses as criteria Sp @Sn=100%.
%      options.b can be any Balu classifier structure (see example).
%      options.show = 1 display results.
%      selec is the indices of the selected features.
%
clt
load datareal                  % f 200 samples x 279 features
op.m = 10;                     % 10 features will be selected
op.show = 1;                   % display results
op.b.name = 'fisher';          % SFS with Fisher
s = Bfs_sfs(f,d,op);           % index of selected features
X = f(:,s);                    % selected features
Xn = fn(s,:)                   % list of feature names
op_lda.p = [];
ds = Bcl_lda(X,d,X,op_lda);    % LDA classifier
p = Bev_performance(d,ds)      % performance with sfs


%% 3.2 FEATURE SELECTION - RANK
% selec = Bfs_rank(X,d,options)
%
%      Feature selection based on command rankfeatures (from MATLAB
%      Bioinformatics Toolbox) that ranks ranks key features by class
%      separability criteria.
%
%      input: X feature matrix
%             options.m number of features to be selected
%             options.criterion can be:
%        'ttest' (default) Absolute value two-sample T-test with pooled
%                          variance estimate
%        'entropy'         Relative entropy, also known as Kullback-Lieber
%                          distance or divergence
%        'brattacharyya'   Minimum attainable classification error or
%                          Chernoff bound
%        'roc'             Area between the empirical receiver operating
%                          characteristic (ROC) curve and the random classifier
%                          slope
%        'wilcoxon'        Absolute value of the u-statistic of a two-sample
%                          unpaired Wilcoxon test, also known as Mann-Whitney
%
%       Notes: 1) 'ttest', 'entropy', and 'brattacharyya' assume normal
%       distributed classes while 'roc' and 'wilcoxon' are nonparametric tests,
%       2) all tests are feature independent.
%
%      output: selec selected features
clt
load datareal
op.m = 50;                     % 10 features will be selected
op.criterion = 'roc';          % ROC criterion will be used
op.show = 1;                   % display results
s = Bfs_rank(f,d,op);          % index of selected features
X = f(:,s);                    % selected features
Xn = fn(s,:)                   % list of feature names
op_knn.k = 3;
ds = Bcl_knn(X,d,X,op_knn);    % KNN classifier
p = Bev_performance(d,ds)      % performance with rank


%% 4.1 FEATURE TRANSFORMATION - PCA
% [Y,lambda,A,Xs,mx] = Bft_pca(X,m)
%
%      Principal component analysis
%      X is the matrix feature.
%      m number of selected components or the energy 0<m<=1 (in this case it
%      will be selected the first d principal components that fulfill the
%      condition sum(lambda(1:d))/sum(lambda) >= energy. See Example 2)
%      Y contains the m principal components of X
%      lambda contains the standard deviation of each principal component
%      A is the transformation matrix from X to Y
%      Xs is the reconstructed X from A and Y.
%
clt
X = double(imread('cameraman.tif')); % 256x256 pixels
[Y,lambda,A,Xs] = Bft_pca(X,30);     % 30 principal components
figure(1);bar(lambda/lambda(1))
figure(2);imshow([X Xs],[])


%% 4.2 FEATURE TRANSFORMATION - PLSR
% [T,U,P,Q,W,B] = Bft_plsr(X,d,options)
%
%      Feature transformation using Partial Least Squares Regression with
%      NIPALS algorithm.
%      X: Input matrix with features
%      d: Vector with ideal classifcation.
%      m: Number of principal components to be selected.
%      T: Loadings of X (m transformed features). Matrix T is Xo*W', where
%         Xo is normalized X and W corresponds to the weights.
%      U: Loadings of Y (m transformed features). Matrix U is T*B.
%      P: Principal components of X
%      Q: Principal components of Y (defined from d, one binary column
%         per class)
%      W: X weights
%      B: A vector of regression coefficients
clt
load datareal
s = [279 235 268 230 175 165 207 160 269 157]; %indices using Example of Bfs_sfs.
X = f(:,s);
T = Bft_plsr(X,d,6);
op.p = [];
ds1 = Bcl_lda(T,d,T,op);
Bev_performance(d,ds1)

% Comparison with first 6 SFS features
ds2 = Bcl_lda(X(:,1:6),d,X(:,1:6),op);
Bev_performance(d,ds2)



%% 5.1 CLASSIFICATION USING LDA
%   ds      = Bcl_lda(X,d,Xt,[])  Training & Testing together
%   options = Bcl_lda(X,d,[])     Training only
%   ds      = Bcl_lda(Xt,options) Testing only
%
%      LDA (linear discriminant analysis) classifier.
%      We assume that the classes have a common covariance matrix
%
%      Design data:
%         X is a matrix with features (columns)
%         d is the ideal classification for X
%         options.p is the prior probability, if p is empty,
%         it will be estimated proportional to the number of samples of each
%         class.
%
%      Test data:
%         Xt is a matrix with features (columns)
%
%      Output:
%         ds is the classification on test data
%         options.dmin contains min(d).
%         options.Cw1 is inv(within-class covariance).
%         options.mc contains the centroids of each class.
%         options.string is a 8 character string that describes the performed
%         classification (in this case 'lda     ').
%
%
%      Example: Training & Test together:
clt
load datagauss             % simulated data (2 classes, 2 features)
Bio_plotfeatures(X,d)      % plot feature space
op.p = [];
ds = Bcl_lda(X,d,Xt,op);   % LDA classifier
p = Bev_performance(ds,dt) % performance on test data

%% 5.2 CLASSIFICATION USING BAYES
%   ds      = Bcl_bayes2(X,d,Xt,options)  Training & Testing together
%   options = Bcl_bayes2(X,d,options)     Training only
%   ds      = Bcl_bayes2(Xt,options)      Testing only
%
%      Bayes classifier for ONLY two features and two classes
%
%      Design data:
%         X is a matrix with features (columns)
%         d is the ideal classification for X
%         options.p is the prior probability, if p is not given, it will be estimated
%           proportional to the number of samples of each class.
%         options.show = 1 displays results.
%
%      Test data:
%         Xt is a matrix with features (columns)
%
%      Output:
%         ds is the classification on test data
%         options.dmin contains min(d).
%         options.string is a 8 character string that describes the performed
%         classification (in this case 'bayes2  ').
%         options.m1,m2,b1,b2,C are parameters of the classifier.

load datagauss              % simulated data (2 classes, 2 features)
Bio_plotfeatures(X,d)       % plot feature space
op.p = [0.25 0.75];         % prior probability
op.show = 1;                % display results
ds = Bcl_bayes2(X,d,Xt,op); % Bayes classifier
p = Bev_performance(ds,dt)  % performance on test data


%% 5.3 CLASSIFICATION USING RANDOM FOREST
%   ds      = Bcl_RandomForest(X,d,Xt,[])  Training & Testing together
%   options = Bcl_RandomForest(X,d,[])     Training only
%   ds      = Bcl_RandomForest(Xt,options) Testing only
%
%      Classifier using Random Forest. This implementation uses command
%      TreeBagger of Statistics and Machine Learning Toolbox of Matworks
%
%      Design data:
%         X is a matrix with features (columns)
%         d is the ideal classification for X
%
%      Test data:
%         Xt is a matrix with features (columns)
%
%      Output:
%         ds is the classification on test data
%         options.NTrees number of trees
%         options.string is a 8 character string that describes the performed
%         classification (in this case 'RForest ').
%
%      Example: Training & Test together:
clt
load datagauss                     % simulated data (2 classes, 2 features)
Bio_plotfeatures(X,d)              % plot feature space
op.NTrees  = 50;                   % Number of trees
ds = Bcl_RandomForest(X,d,Xt,op);  % RandomForest classifier
p = Bev_performance(ds,dt)         % performance on test data


%% 5.4 CLASSIFICATION USING SVM
%   ds      = Bcl_svm(X,d,Xt,options)  Training & Testing together
%   options = Bcl_svm(X,d,options)     Training only
%   ds      = Bcl_svm(Xt,options)      Testing only
%
%   Toolbox: Balu
%      Support Vector Machine approach using the LIBSVM(*).
%
%      Design data:
%         X is a matrix with features (columns)
%         d is the ideal classification for X
%
%         options.kernel is a string that defines the options of LIBSVM as
%         follows:
%
%            -s svm_type : set type of SVM (default 0)
%                    0 -- C-SVC
%                    1 -- nu-SVC
%                    2 -- one-class SVM
%                    3 -- epsilon-SVR
%                    4 -- nu-SVR
%            -t kernel_type : set type of kernel function (default 2)
%                    0 -- linear: u'*v
%                    1 -- polynomial: (gamma*u'*v + coef0)^degree
%                    2 -- radial basis function: exp(-gamma*|u-v|^2)
%                    3 -- sigmoid: tanh(gamma*u'*v + coef0)
%            -d degree : set degree in kernel function (default 3)
%            -g gamma : set gamma in kernel function (default 1/num_features)
%            -r coef0 : set coef0 in kernel function (default 0)
%            -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
%            -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
%            -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
%            -m cachesize : set cache memory size in MB (default 100)
%            -e epsilon : set tolerance of termination criterion (default 0.001)
%            -h shrinking: whether to use the shrinking heuristics, 0 or 1 (default 1)
%            -b probability_estimates: whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
%            -wi weight: set the parameter C of class i to weight*C, for C-SVC (default 1)%
%
%      Test data:
%         Xt is a matrix with features (columns)
%
%      Output:
%         ds is the classification on test data
%         options.svmStruct contains information about the trained classifier
%         (from function svmtrain of Bioinformatics Toolbox).
%         options.string is a 8 character string that describes the performed
%         classification (e.g., 'Bsvm,4  ' means rbf-SVM).
%
%      Example: Training & Test together:
clt
load datagauss             % simulated data (2 classes, 2 features)
Bio_plotfeatures(X,d)      % plot feature space
op.kernel = '-t 2';
ds = Bcl_libsvm(X,d,Xt,op);% rbf-SVM classifier
p = Bev_performance(ds,dt) % performance on test data

%      http://www.csie.ntu.edu.tw/~cjlin/libsvm
%
%      WARNING: When installing libsvm you must rename the compiled files
%      svmtrain and svmpredict (in libsvmxxx/matlab folder) as
%      libsvmtrain and libsvmpredict. Additionally, libsvmxxx/matlab folder
%      must be included in the matlab path.
%

%% 5.5: CLASSIFICATION - SEVERAL CLASSIFIERS
% ds      = Bcl_structure(X,d,Xt,options)  Training & Testing together
%   options = Bcl_structure(X,d,options)     Training only
%   ds      = Bcl_structure(Xt,options)      Testing only
%
%      Classification using Balu classifier(s) defined in structure b.
%
%      Design data:
%         X is a matrix with features (columns)
%         d is the ideal classification for X
%         options is a Balu classifier structure b with
%            b.name      = Balu classifier's name
%            b.options   = options of the classifier
%
%         b can define one or more classifiers (see example).
%
%      Test data:
%         Xt is a matrix with features (columns)
%
%      Output:
%         ds is the classification on test data (one column per classifier)
%
%      Example: Training & Test together:
clt
load datagauss                                                        % simulated data (2 classes, 2 features)
b(1).name = 'knn';   b(1).options.k = 5;                              % KNN with 5 neighbors
b(2).name = 'knn';   b(2).options.k = 7;                              % KNN with 7 neighbors
b(3).name = 'knn';   b(3).options.k = 9;                              % KNN with 9 neighbors
b(4).name = 'lda';   b(4).options.p = [];                             % LDA
b(5).name = 'qda';   b(5).options.p = [];                             % QDA
b(6).name = 'nnglm'; b(6).options.method = 3; b(6).options.iter = 10; % Nueral network
b(7).name = 'svm';   b(7).options.kernel = 4;                         % rbf-SVM
b(8).name = 'maha';  b(8).options = [];                               % Euclidean distance
b(9).name = 'dmin';  b(9).options = [];                               % Mahalanobis distance
op = b;
ds = Bcl_structure(X,d,Xt,op);                                        % ds has 9 columns
p = Bev_performance(ds,dt)

for i=1:length(b)
    fprintf('%10s: %7.4f\n',b(i).name,p(i));
end


%% 5.6 CLASSIFICATION - ENSEMBLE
%   ds      = Bcl_ensemble(X,d,Xt,options)  Training & Testing together
%   options = Bcl_ensemble(X,d,options)     Training only
%   ds      = Bcl_ensemble(Xt,options)      Testing only
%
%      Design and test an ensemble of n classifiers.
%
%      Design data:
%         X is a matrix with features (columns)
%         d is the ideal classification for X
%         tensemble is the type of the ensemble:
%            'vote'  for majority vote
%            'min'   for minimum
%            'max'   for maximum
%            'best'  takes the best individual classifier reclassifying X
%            'all'   for all classifications (one classification per column in ds)
%            'stack' uses a stacked generalization acoording to Polikar-Fig 11.
%                    in this case param is the classifier defined as b
%                    structure.
%
%         b is the following structure for classifier i (i=1,...,n)
%            b(i).name = classifier's name
%            b(i).p1   = parameter 1 of the classifier if any
%            b(i).p2   = parameter 2 of the classifier if any
%            b(i).p3   = parameter 3 of the classifier if any
%
%      Test data:
%         Xt is a matrix with features (columns)
%
%      Output:
%         ds is the classification on test data
%
%      show = 1 displays results (default show = 0).
%
%      Example: Training & Test together:
clt
load datagauss
b(1).name = 'knn';   b(1).options.k = 5;                              % KNN with 5 neighbors
b(2).name = 'knn';   b(2).options.k = 7;                              % KNN with 7 neighbors
b(3).name = 'knn';   b(3).options.k = 9;                              % KNN with 9 neighbors
b(4).name = 'lda';   b(4).options.p = [];                             % LDA
b(5).name = 'qda';   b(5).options.p = [];                             % QDA
b(6).name = 'nnglm'; b(6).options.method = 3; b(6).options.iter = 10; % Nueral network
b(7).name = 'svm';   b(7).options.kernel = 4;                         % rbf-SVM
b(8).name = 'maha';  b(8).options = [];                               % Euclidean distance
b(9).name = 'dmin';  b(9).options = [];                               % Mahalanobis distance
op.b   = b;
op.tensemble = 'vote';
op.show = 1;
op.bs   = [];
ds1 = Bcl_ensemble(X,d,Xt,op);               % majority vote
p1  = Bev_performance(ds1,dt)                % performance on test data

op.tensemble = 'max';
op.bs   = [7 0.8 0.7];
ds2 = Bcl_ensemble(X,d,Xt,op);               % DCS
p2  = Bev_performance(ds2,dt)                % performance on test data


%% 5.7: CLASSIFICATION SELECTION USING GUI
%
%
clt
cd orl
load Bfx_results
d = double(Bds_labels(10*ones(40,1)));
cd ..
save fx_orl f fn d
disp('In the GUI please select 40 features, SFS, LDA, fx_orl.mat file and GO!');
enterpause
Bcl_gui2


%% 6.1 MULTIPLE VIEWS - FUNDAMENTAL MATRIX
% F = Bmv_fundamental(A,B,method)
%
%
%   Fundamental matrix from projection matrices.
%
%   F = fundamental(A,B,method) returns the 3x3 fundamental matrix from
%      3x4 projection matrices A and B according to the following
%      methods:
%
%      method = 'tensor' : uses bifocal tensors with canonic matrix
%
%      method = 'tensor0': uses bifocal tensors without canonic matrix
%
%      method = 'pseudo' : uses pseudoinverse matrix of A (deault)
%
%      If no method is given, 'pseudo' will be assumed as default.
%
%      Both methods can be found in:
%
%      R. Hartley and A. Zisserman. Multiple View Geometry in Computer
%      Vision. Cambridge University Press, 2000.
%
%      Example:
clt
A = rand(3,4);                      % Projection matrix for view 1
B = rand(3,4);                      % Projection matrix for view 2
M = [1 2 3 1]';                     % 3D point (X=1,Y=2,Z=3)
m1 = A*M;m1 = m1/m1(3);             % projection point in view 1
m2 = B*M;m2 = m2/m2(3);             % projection point in view 2
F = Bmv_fundamental(A,B,'tensor');  % Fundamental matrix using tensors
m2'*F*m1                            % epipolar constraint must be zero

%% 6.2 MULTIPLE VIEWS - TRIFOCAL TENSOR
% F = Bmv_trifocal(A,B,C)
%
%      Trifocal tensors.
%
%      T = trifocal(A,B,method) returns the trifocal tensors from
%      3x4 projection matrices A, B and C according of thre
%      views. T is a 3x3x3 array
%
clt
A = rand(3,4);           % Projection matrix for view 1
B = rand(3,4);           % Projection matrix for view 2
C = rand(3,4);           % Projection matrix for view 3
T = Bmv_trifocal(A,B,C)  % Trifocal tensors


%% 6.3 MULTIPLE VIEWS - ESTIMATION OF FUNDAMENTAL MATRIX USING SIFT
% F = Bmv_fundamentalSIFT(I1,I2)
%
%
%      Estimation of Fundamental Matrix from two images using SIFT points.
%
%      This function requires VLFeat Toolbox from (www.vlfeat.org).
%
%      I1 and I2 are the stereo images.
%      F is the Fundamental matrix.
%
clt
I1 = imread('testimg5.jpg');            % Image 1
figure(1);imshow(I1); hold on
I2 = imread('testimg6.jpg');            % Image 2
figure(2);imshow(I2); hold on
F  = Bmv_fundamentalSIFT(I1,I2);        % F matrix estimation
while(1)
    figure(1);
    disp('click a point in Figure 1...') % click
    p = vl_click; m1 = [p(2) p(1) 1]';
    plot(p(1),p(2),'g+')
    figure(2)
    Bmv_epiplot(F,m1)
end


%% 6.4 MULTIPLE VIEWS - HOMOGRAPHY
%   [Ibs,H] = Bmv_homographySIFT(Ia,Ib,show)
%
%
%      Homography between images Ia and Ib using RANSAC of SIFT points
%      Ibs is transformed image Ib.
%      size(Ia) = size(Ibs)
%      If ma and mb are the homogeneus coordinates of points in images Ia and Ib:
%         ma = [xa ya 1]'; mb = [xb yb 1]';
%         then  mb = H*ma.
%      show = 1 displays results.
%
clt
Ia = rgb2gray(imread('testimg4.jpg'));           % Image Ia
H  = [1.9 0.01 0.01;0.002 1.95 0.001;0.001 0 1]  % Original H
Ib = Bmv_projective2D(Ia,H,[1000 750],1);        % Image Ib
[Ibs,Hs] = Bmv_homographySIFT(Ia,Ib,1);          % Homography
Hs = Hs/Hs(3,3)                                  % estimated matrix H

%% 6.5 MULTIPLE VIEWS - PROJECTIVE TRANSFORMATION
%  J = Bmv_projective2D(I,H,SJ,show)
%
%
%      2D proyective transformation.
%
%      J = projective2D(I,H,SJ,show) returns a new image J that is computed
%      from the 2D projective transformation H of I.
%
%      SJ is [NJ MJ] the size of the transformed image J. The
%      default of SJ is [NJ,MJ] = size(I).
%
%      show = 1 displays images I and J.
%
%      The coordinates of I are (xp,yp) (from x',y'), and the
%      coordinates of J are (x,y). The relation between both
%      coordinate systems is: mp = H*m, where mp = [xp yp 1]'
%      and m = [x y 1]'.
%
clt
I = rgb2gray(imread('testimg4.jpg'));  % original image
s = 7.5; t = pi/4;                     % scale and orientation
H = [s*cos(t) -s*sin(t)  -500
    s*sin(t)  s*cos(t) -1750
    0         0        1];         % similarity matrix
J = Bmv_projective2D(I,H,[600 250],1); % projective transformation


%% 7.1 CLUSTERING - KMEANS
%  function [ds,C] = Bct_kmeans(X,k,show)
%
%
%      k-means clustering
%      X matrix of samples
%      k number of clusters
%      ds assigned class number
%      C centroids of each cluster
%      show = 1 display intermediate results
%      show = 0 uses kmeans algortithm of vlfeat (if installed else algorithm
%      of matlab).
%
clt
[X,d] = Bds_gaussgen([10 1;1 10;15 15],4*ones(3,3),100*ones(3,1));
figure(1)
Bio_plotfeatures(X,d);
ds = Bct_kmeans(X,3);
figure(2)
Bio_plotfeatures(X,ds);

%% 7.2 CLUSTERING - MEANSHIFT
%   function [ds,Z] = Bct_meanshift(X, sigma)
%
%      Meanshift clustering
%      X matrix of samples
%      sigma: standard deviation of the Gaussian Parzen window
%      ds assigned class number.
%      Z are the reduced coordinates.
%
clt
[X,d] = Bds_gaussgen([10 1;1 10;15 15],4*ones(3,3),100*ones(3,1));
figure(1)
Bio_plotfeatures(X,d);
ds = Bct_meanshift(X,3);
figure(2)
Bio_plotfeatures(X,ds);


%% 7.3 CLUSTERING - NEIGHBOR
%   function [ds,Xc] = Bct_neighbor(X,th)
%
%      Neigbor clustering: iterative method, a sample will be added to a
%      cluster if its distance to the mass center of the cluster is less than
%      th, else it will be create a new cluster.
%
%      X matrix of samples
%      th minimal distance
%      ds assigned class number
%      Xc mass center of each cluster
%
clt
[X,d] = Bds_gaussgen([10 1;1 10],1*ones(2,2),100*ones(2,1));
figure(1)
Bio_plotfeatures(X,d);
ds = Bct_neighbor(X,4);
figure(2)
Bio_plotfeatures(X,ds);

%% 8.1 EVALUATION - HOLDOUT
%   [T,p] = Bev_holdout(X,d,options)
%
%      Holdout evaluation of a classifier.
%
%      X is a matrix with features (columns)
%      d is the ideal classification for X
%
%      options.b is a Balu classifier or several classifiers (see example)
%      options.s is the portion of data used for training, e.g. s=0.75.
%      options.strat = 1 means the portion is stratified.
%
%      Confusion matrix is given file by file for each group in matrix T.
%      The performance in each group is given vector per.
%
%      Example for one classifier:
clt
load datagauss                      % simulated data (2 classes, 2 features)
Bio_plotfeatures(X,d)                  % plot feature space
b.name = 'knn'; b.options.k = 5;   % knn with 5 neighbors
op.b = b; op.strat = 1;op.s = 0.75; % stratify with 75% train
[T,p] = Bev_holdout(X,d,op)         % holdout
enterpause
%     Example for more classifiers:
load datagauss                                                        % simulated data (2 classes, 2 features)
b(1).name = 'knn';   b(1).options.k = 5;                              % KNN with 5 neighbors
b(2).name = 'knn';   b(2).options.k = 7;                              % KNN with 7 neighbors
b(3).name = 'knn';   b(3).options.k = 9;                              % KNN with 9 neighbors
b(4).name = 'lda';   b(4).options.p = [];                             % LDA
b(5).name = 'qda';   b(5).options.p = [];                             % QDA
b(6).name = 'nnglm'; b(6).options.method = 3; b(6).options.iter = 10; % Nueral network
b(7).name = 'svm';   b(7).options.kernel = 4;                         % rbf-SVM
b(8).name = 'maha';  b(8).options = [];                               % Euclidean distance
b(9).name = 'dmin';  b(9).options = [];                               % Mahalanobis distance
op.b = b; op.strat = 1; op.s = 0.75; op.show=1;                       % stratify with 75% train
[T,p] = Bev_holdout(X,d,op);                                          % holdout

%% 8.2 EVALUATION - CROSS-VALIDATION
%  [p,ci] = Bev_crossval(X,d,options)
%
%
%      Cross-validation evaluation of a classifier.
%
%      v-fold Cross Validation in v groups of samples X and classification d
%      according to given method. If v is equal to the number of samples,
%      i.e., v = size(X,1), this method works as the original
%      cross-validation, where training will be in X without sample i and
%      testing in sample i. ci is the confidence interval.
%
%      X is a matrix with features (columns)
%      d is the ideal classification for X
%
%      options.b is a Balu classifier or several classifiers (see example)
%      options.v is the number of groups (folders) of the cross-validations
%      options.c is the probability of the confidence intervale.
%      options.strat = 1 means the portions are stratified.
%      options.show = 1 displays results.
%
%      Example for one classifier without stratified data:
clt
load datagauss                                                        % simulated data (2 classes, 2 features)
Bio_plotfeatures(X,d)                                                    % plot feature space
b.name = 'knn'; b.options.k = 5;                                      % knn with 5 neighbors
op.strat = 0; op.b = b; op.v = 10; op.c = 0.90; op.show = 0;     	  % 10 groups cross-validation for 90%
[p,ci] = Bev_crossval(X,d,op)                                         % cross valitadion
enterpause

%     Example for more classifiers with stratified data:
load datagauss                                                        % simulated data (2 classes, 2 features)
b(1).name = 'knn';   b(1).options.k = 5;                              % KNN with 5 neighbors
b(2).name = 'knn';   b(2).options.k = 7;                              % KNN with 7 neighbors
b(3).name = 'knn';   b(3).options.k = 9;                              % KNN with 9 neighbors
b(4).name = 'lda';   b(4).options.p = [];                             % LDA
b(5).name = 'qda';   b(5).options.p = [];                             % QDA
b(6).name = 'nnglm'; b(6).options.method = 3; b(6).options.iter = 10; % Nueral network
b(7).name = 'libsvm';b(7).options.kernel = '-t 2';                    % rbf-SVM
b(8).name = 'dmin';  b(8).options = [];                               % Euclidean distance
b(9).name = 'maha';  b(9).options = [];                               % Mahalanobis distance
op.strat=1; op.b = b; op.v = 10; op.show = 1; op.c = 0.95;        	  % 10 groups cross-validation
[p,ci] = Bev_crossval(X,d,op);                                        % cross valitadion

%% 8.3 EVALUATION - ROC
%   [Az,Sn,Sp1,t] = Bev_roc(z,d,show)
%
%      ROC Analysis for feature z with classification c. show = 1 indicates
%      that the ROC curve will be displayed.
%      Az is the area under the ROC curve
%      Sn and Sp1 are the coordinates of Sensitibity and 1-Specificity of
%      optimal point of ROC curve (nearest point to (0,1).
%      t is the threshold for this optimal point.
%
clt
load datagauss          % simulated data (2 classes, 2 features)
Bev_roc(X(:,1),d,1)     % ROC plot of first feature of X

Az = Bev_roc(X(:,2),d); % Az of ROC for second feature of X

%% 9.1 INPUT & OUTPUT - DECISION LINE
%
%   Bio_decisionline(X,d,op)
%
%
%      Diaplay a 2D feature space and decision line.
%
%      X: Sample data
%      d: classification of samples
%      op: output of a trained classifier.
%
clt
load datagauss                             % simulated data (2 classes, 2 features)
Xn = ['\beta_1';'\beta_2'];
b(1).name = 'knn';   b(1).options.k = 5;                              % KNN with 5 neighbors
b(2).name = 'lda';   b(2).options.p = [];                             % LDA
b(3).name = 'svm';   b(3).options.kernel = 4;                         % rbf-SVM
op = b;
op = Bcl_structure(X,d,op);
close all
Bio_decisionline(X,d,Xn,op);

%% 9.2 INPUT & OUTPUT - LABEL REGIONS
%   [d,D] = Bio_labelregion(I,L,c)
%
%      User interface to label regions of an image.
%
%      I is the original image (color or grayvalue).
%      L is a labeled image that indicates the segmented regions of I.
%      c is the maximal number of classes.
%      d(i) will be the class number of region i.
%      D is a binary image with the corresponding labels.
%
%   Example:
I = imread('rice.png');
I = I(1:70,1:70);                              % input image
[L,m] = Bim_segmowgli(I,ones(size(I)),40,1.5); % segmented image
[d,D] = Bio_labelregion(I,L,3);


%% 9.3 INPUT & OUTPUT -
%
%   Bio_plotfeatures(X,d,Xn)
%
%   Toolbox: Balu
%      Plot features X acording classification d. If the feature names are
%      given in Xn then they will labeled in each axis.
%
%      For only one feature, histograms are ploted.
%      For two (or three) features, plots in 2D (or 3D) are given.
%      For m>3 features, m x m 2D plots are given (feature i vs. feature j)
%
%   Example 1: 1D & 2D
clt
load datagauss                    % simulated data (2 classes, 2 features)
figure(1)
Bio_plotfeatures(X(:,1),d,'x1')   % histogram of feature 1
figure(2)
Bio_plotfeatures(X,d)             % plot feature space in 2D (2 features)
enterpause
%  Example 2: 3D
load datareal                     % real data
X = f(:,[221 175 235]);           % only three features are choosen
figure
Bio_plotfeatures(X,d)             % plot feature space in 3D (3 features)
enterpause
%  Example 3: 5D (using feature selection)
load datareal                     % real data
op.m = 5;                         % 5 features will be selected
op.s = 0.75;                      % only 75% of sample will be used
op.show = 0;                      % display results
op.b.name = 'fisher';             % definition SFS with Fisher
s = Bfs_balu(f,d,op);             % feature selection
figure
Bio_plotfeatures(f(:,s),d)        % plot feature space for 5 features



%% 9.4 INPUT & OUTPUT - SHOW CONFUSION MATRIX
%  function Bio_showconfusion(C)
%
%   Toolbox: Balu
%      Show confusion matrix C in a color 2D representation.
%
%
%   Example:
%
clt
% Simulation of a 10x10 confussion matrix
C = 1.3*rand(10,10)+2*eye(10);C = C/max(C(:));

Bio_showconfusion(C)

%% 9.5 INPUT & OUTPUT - LATEX TABLE
%  function Bio_latextable(row_names,col_names,fmt,T)
%
%      Code for a latex table.
%
%      row_names is a cell with the names of the rows.
%      col_names is a cell with the names of the columns.
%      fmt is a cell with the format of each column.
%      T is the table.
%
clt
col_names = {'cols','col 1','col 2','col 3','col4'};
row_names = {'row 1','row 2','row 3'};
fmt = {'%5.2f','%6.4f','%3.1f','%7.4f'};
T = rand(3,4);
Bio_latextable(row_names,col_names,fmt,T)


