% IMAGE FILES
images = Bim_build('--enter directory here--','jpg');

% SUPERVISION
labels = Bio_labelimage(images);

% FEATURES
options.fx.b = {'fourierdes','hugeo','flusser','haralick','lbp','gabor'};
options.fx.channels = {'gray','red','green','blue','L*','a*','b*'};
options.fx.segmentation = 'Bim_segbalu';

% FEATURE SELECTION ALGORITHMS
options.fs = {'sfs-lda','sfs-fisher','sfs-knn5','rank-roc','rank-ttest'};               

% CLASSIFIERS
options.cl = {'dmin','maha','lda','qda','knn5','nnglm1','svm1'};

% FEATURE EXTRACTION
[features,names] = Bfx_files(images,options.fx);

% FEATURE AND CLASSIFIER SELECTION 
options.Xn  = names;
options.m   = 20;
[cs,fs] = Bcl_balu(features,labels,options);