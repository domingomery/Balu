% DELETE ALL
clt

% IMAGE FILES DEFINTION
bim = Bim_build('/Users/domingomery/Dropbox/Public/ToolboxBalu/ImagesAndData/fishbones/','png');
%bim = Bim_build('/Users/domingomery/Dropbox/Public/ToolboxBalu/ImagesAndData/chips/','jpg',0.1);
%bim = Bim_build('/Users/domingomery/Matlab/papas60/','JPG',0.1);


% FEATURE EXTRACTION DEFINTION
fx_options.b             = {'lbp'};
fx_options.channels      = {'gray'};
%fx_options.segmentation  = 'Bim_segbalu';

% FEATURE SELECTION DEFINITION
fs_options = {'sfs-lda','sfs-fisher','rank-roc','rank-ttest'};               

% CLASSIFIERS DEFINITIONS
cl_options = {'dmin','maha','lda','qda','svm1'};

% FEATURE EXTRACTION
[X,Xn,S] = Bfx_files(bim,fx_options);

% SUPERVISION
d = [ones(100,1);2*ones(100,1)];
%d = Bio_labelimages(bim);

% FEATURE AND CLASSIFIER SELECTION 
options.Xn  = Xn;
options.m   = 20;
options.fs  = fs_options;
options.cl  = cl_options;

[bcs,selec] = Bcl_balu(X,d,options);


