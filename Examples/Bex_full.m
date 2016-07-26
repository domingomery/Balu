% DELETE ALL
clt

% IMAGE FILES DEFINTION
 bim = Bim_build('/Users/domingomery/Dropbox/Public/ToolboxBalu/ImagesAndData/fishbones/','png');
% bim = Bim_build('/Users/domingomery/Dropbox/Public/ToolboxBalu/ImagesAndData/chips/','jpg',0.1);
%bim = Bim_build('/Users/domingomery/Matlab/papas60/','JPG',0.1);


% GEOMETRIC FEATURE EXTRACTION DEFINTION
opg.b            = Bfx_build();

%opg.b            = Bfx_build('allgeo');
%opg.segmentation = 'Bim_segbalu';


% INTENSITY FEATURE EXTRACTION DEFINTIONS
opi.b            = Bfx_build('lbp');
opi.channels      = 'g';
% opi.colstr       = 'gRGBHSV';
% opi.segmentation = 'Bim_segbalu';

% FEATURE SELECTION DEFINITION
bfs = Bfs_build('sfs-lda','sfs-qda','sfs-nnglm2','rank-roc');               

% CLASSIFIERS DEFINITIONS
bcl = Bcl_build('svm1','svm2','adaboost','boosting');

% FEATURE EXTRACTION
[X,Xn,S] = Bfx_files(bim,opg,opi);


% SUPERVISION
d = [ones(100,1);2*ones(100,1)];
%d = Bio_labelimages(bim);

% FEATURE AND CLASSIFIER SELECTION 
options.Xn  = Xn;
options.m   = 20;
options.v   = 10;

[bcs,selec] = Bcl_balu(X,d,bcl,bfs,options);



