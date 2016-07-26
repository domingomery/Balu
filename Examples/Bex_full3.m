% % DELETE ALL
% clt
% 
% % IMAGE FILES DEFINTION
% bim = Bim_build('/Users/domingomery/Dropbox/Public/ToolboxBalu/ImagesAndData/fishbones/','png');
% %bim = Bim_build('/Users/domingomery/Dropbox/Public/ToolboxBalu/ImagesAndData/chips/','jpg',0.1);
% %bim = Bim_build('/Users/domingomery/Matlab/papas60/','JPG',0.1);
% 
% 
% % FEATURE EXTRACTION DEFINTION
% options.fx.b             = {'lbp'};
% options.fx.channels      = {'gray'};
% %fx_options.segmentation  = 'Bim_segbalu';
% 
% % FEATURE SELECTION DEFINITION
% options.fs = {'sfs-lda','sfs-fisher','rank-roc','rank-ttest'};               
% 
% % CLASSIFIERS DEFINITIONS
% options.cl = {'dmin','maha','lda','qda','svm1'};
% 
% % FEATURE EXTRACTION
% [X,Xn,S] = Bfx_files(bim,options.fx);
% 
% % SUPERVISION
% d = [ones(100,1);2*ones(100,1)];
% %d = Bio_labelimages(bim);
% 
% % FEATURE AND CLASSIFIER SELECTION 
% options.Xn  = Xn;
% options.m   = 20;
% [bcs,selec] = Bcl_balu(X,d,options);
% 
% 
% 
% 
clt
% IMAGE FILES
%images = Bim_build('/Users/domingomery/Dropbox/Public/ToolboxBalu/ImagesAndData/chips/','jpg',0.1);
images = Bim_build('/Users/domingomery/Dropbox/Public/ToolboxBalu/ImagesAndData/fishbones/','png');

% SUPERVISION
% labels = Bio_labelimages(images);
labels = [ones(100,1);2*ones(100,1)];;
% FEATURES
%options.fx.b = {'fourierdes','hugeo','flusser','haralick','lbp'};
%options.fx.channels = {'gray','red','green','blue'};
%options.fx.segmentation = 'Bim_segbalu';

options.fx.b = {'lbps'};
options.fx.channels = {'gray'};




% FEATURE SELECTION ALGORITHMS
options.fs = {'sfs-lda','sfs-fisher','sfs-knn5','rank-roc','rank-ttest'};               

% CLASSIFIERS
options.cl = {'dmin','maha','lda','qda','knn5','nnglm1','svm1'};

% FEATURE EXTRACTION
[features,names] = Bfx_files(images,options.fx);

% FEATURE AND CLASSIFIER SELECTION 
options.Xn  = names;
options.m   = 20;
[bcs,selec] = Bcl_balu(features,labels,options);