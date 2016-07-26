% Contents.m
%
% Toolbox: Balu
%    Balu Help
%
%    Display help information of Balu Toolbox Matlab
%
% D.Mery 2014
% http://dmery.ing.puc.cl
%
% Image processing (im)
%    Bim_build             + Build image structure for Bio_loadimg
%    Bim_color2bwreg       + Convert a region to b&w
%    Bim_colorenhancement  + Enhancement of a color image
%    Bim_colorconv         + Color space conversion
%    Bim_cssalient         + Saliency map
%    Bim_d1                + First derivative (gradient)
%    Bim_d2                + Second derivative (laplacian)
%    Bim_deconvolution     + Inverse filtering
%    Bim_equalization      + Equalization forcing uniform histogram
%    Bim_fconversion       + image format conversion
%    Bim_hsi2rgb           + Conversion HSI to RGB
%    Bim_inthist           + Integral histogram
%    Bim_inthistread       + Integral histogram of a part of an image
%    Bim_labparam          + Estimation of L*a*b* parameters
%    Bim_lin               + Lineal enhancement
%    Bim_LUT               + Grayvalue transformation using lookup table
%    Bim_maxmin            + Normalization of a image between 0 and 1
%    Bim_morphoreg         + Morphological operations
%    Bim_performance       + Precision - recall using real and ideal segmentation  
%    Bim_regiongrow        + Interactive segmentation
%    Bim_resminio          + Restoration using MINIO criterium
%    Bim_rgb2hcm           + Conversion RGB to high contrast image
%    Bim_rgb2hsi           + Conversion RGB to HSI
%    Bim_rgb2lab           + Conversion RGB to L*a*b* using calibration
%    Bim_rgb2lab0          + Conversion RGB to L*a*b* using formulas
%    Bim_rgb2pca           + Conversion RGB to PCA
%    Bim_sat               + Saturation of an image
%    Bim_segbalu           + Segmentation of an object using hci transform
%    Bim_segblops          + Segmentation of blops
%    Bim_segdefects        + Segmentation of defects
%    Bim_segkmeans         + Kmeans Segmentation 
%    Bim_segmaxfisher      + Color segmentacion using Fisher discriminat
%    Bim_segmaxvar         + Color segmentacion maximizing the variance
%    Bim_segmowgli         + Segmentation of an object's region
%    Bim_segmser           + Maximal Stable Energy Region Segmentation
%    Bim_segotsu           + Otsu Segmentation 
%    Bim_segpca            + Segmentation of an object using pca
%    Bim_segsliwin         + Segmentation using sliding windows 
%    
% Feature extraction (fx)
%    Bfx_all               + (int) grayvalues of all pixels
%    Bfx_basicgeo          + (geo) Standard geometric features
%    Bfx_basicint          + (int) Standard intesnity features
%    Bfx_build             + Construction of a Balu Feature Extraction structure
%    Bfx_clp               + (int) Crossing line profiles
%    Bfx_contrast          + (int) Contrast features
%    Bfx_dct               + (int) DCT features
%    Bfx_files             + (int & geo) feature extraction from a set of images
%    Bfx_fitellipse        + (geo) Ellipse features
%    Bfx_flusser           + (geo) Flusser moments
%    Bfx_fourier           + (int) Fourier features
%    Bfx_fourierdes        + (geo) Fourier descriptors
%    Bfx_gabor             + (int) Gabor texture features
%    Bfx_geo               + (geo) Geometric features
%    Bfx_gui               + Graphic user interface for feature extraction
%    Bfx_gupta             + (geo) Gupta moments
%    Bfx_haralick          + (int) Haralick texture features
%    Bfx_hog               + (int) Histogramf of Oriented Gradients
%    Bfx_hugeo             + (geo) Hu moments
%    Bfx_huint             + (int) Hu moments with intensity
%    Bfx_int               + (int) Intensity features
%    Bfx_lbp               + (int) Local binary pattern texture features (6)
%    Bfx_lbpcontrast       + (int) Contrast using LBP features
%    Bfx_lbhog             + (int) LBP and HOG features
%    Bfx_moments           + (geo) Statistical and central moments
%    Bfx_onesift           + (int) Extract one sift descriptor of a region
%    Bfx_phog              + (int) Pyramid Histogram of Oriented Gradients (7)
%    Bfx_randomsliwin      + (int) fx from random sliding windows
%    Bfx_vlhog             + (int) HOG using VLfeat Toolbox
% 
% Feature transformation (ft)
%    Bft_lseft             + Transformation using LSEF algorithm
%    Bft_norm              + Normalization of features
%    Bft_pca               + Principal Components Analysis
%    Bft_pcr               + Principal Components Regression
%    Bft_plsr              + Partial Least Squares Regression (8)
%    Bft_uninorm           + Norm of each column is one
%    Bft_vq                + Vector quantization
% 
% Feature selection (fs)
%    Bfs_all               + Dummy selection: it selects all features
%    Bfs_balu              + Normalize, clean and select features
%    Bfs_bb                + Branch and Bound feature selection
%    Bfs_build             + Construction of a Balu Feature Selection structure
%    Bfs_clean             + Delete high correlated and constant features
%    Bfs_exsearch          + Feature selection using exhaustive serach
%    Bfs_fosmod            + Feature Selection using FOS-MOD algorithm
%    Bfs_lsef              + Feature Selection using LSE-forward algorithm
%    Bfs_mRMR              + Feature selection after Peng
%    Bfs_nobackground      + Delete contrast features
%    Bfs_noposition        + Delete position features
%    Bfs_norotation        + Delete rotation variant features
%    Bfs_random            + Select best features from random subsets
%    Bfs_rank              + Feature selection using command rankfeatures
%    Bfs_ransac            + Select best features from random samples
%    Bfs_sfs               + Sequential forward feature selection
%    Bfs_sfscorr           + SFS for regression problems  
%
% Input/Output (io)
%    Bio_copyfiles         + Copy files of a directory wit a new name
%    Bio_decisionline      + Line decision for 2 features problem
%    Bio_drawellipse       + Draw an ellipse
%    Bio_edgeview          + Edge of an object
%    Bio_findex            + Index number of a feature for a given string
%    Bio_fmtconv           + Image format conversion
%    Bio_imgshow           + Display an image of a sequence
%    Bio_labelimage        + User interface to label a set of images
%    Bio_labelregion       + User interface to label regions of an image
%    Bio_latextable        + Build LateX table
%    Bio_loadimg           + Load an image of a sequence
%    Bio_maillist          + Send emails to a mail-list
%    Bio_plotfeatures      + Feature space
%    Bio_plotroc           + Plot ROC 
%    Bio_plotsquare        + Square onto an image
%    Bio_printfeatures     + Feature values
%    Bio_regshow           + Display a region in a color image
%    Bio_segshow           + Display image and segmentation with Bim_segbalu
%    Bio_sendmail          + Send e-mail
%    Bio_showconfusion     + Show confusion matrix using colormap
%    Bio_statusbar         + Display progress bar (10)
%
% Data selection & generation (ds)
%    Bds_Adasyn            + Oversampling method Adasyn
%    Bds_bootstrap         + Bootstrap sample
%    Bds_BorderSMOTE       + Oversampling method Border SMOTE
%    Bds_CNNRule           + Undersampling method CNN rule
%    Bds_gaussgen          + Generation of Gaussian Data
%    Bds_ixstratify        + Data Stratification (without replacement)
%    Bds_labels            + Generation of (supervised) label vector
%    Bds_NCLRule           + Undersampling method NCL
%    Bds_nostratify        + Data Sampling without Stratification
%    Bds_OSS               + Undersampling method OSS
%    Bds_ROS               + Oversampling method ROS
%    Bds_RUS               + Undersampling method RUS
%    Bds_stratify          + Data Stratification
%    Bds_smote             + Oversampling method SMOTE
%    Bds_TomekLinks        + Undersampling method TomeLinks
% 
% Classification (cl)
%    Bcl_adaboost          + AdaBoost.M2
%    Bcl_AdaboostM1        + AdaBoost.M1
%    Bcl_bagging           + Bagging 
%    BclBalanceCasacade    + Balance Cascade
%    Bcl_balu              + Exhaustive Feature & Classifier Selection
%    Bcl_bayes2            + Bayes for two features
%    Bcl_boosting          + Boosting
%    Bcl_boostVJ           + Viola-Jones boosting
%    Bcl_build             + Construction of a Balu classifier structure
%    Bcl_construct         + Inner procedure (not a classifier function)
%    Bcl_dcs               + Dynamic classifier selection
%    Bcl_det21             + Linear Detector for 2 features
%    Bcl_det22             + Quardatic Detector for 2 features
%    Bcl_dmin              + Euclidean minimal distance
%    Bcl_EasyEnsemble      + Easy Ensemble
%    Bcl_ensemble          + Ensemble of n classifiers
%    Bcl_exe               + Execute a Balu classifier
%    Bcl_gui               + Graphic user interface for model selection
%    Bcl_knn               + k-nearest neighbors (1)
%    Bcl_lda               + Linear discriminant analysis
%    Bcl_libsvm            + SVM using library LIBSVM (9)
%    Bcl_maha              + Mahalanobis minimal distance 
%    Bcl_nbnnxi            + Naive Bayes Nearest Neighbor for histograms
%    Bcl_nnglm             + Neural Network (5)
%    Bcl_outscore          + Function for include score into class vector
%    Bcl_pegasos           + Pegasos Support Vector Machine (1)
%    Bcl_pnn               + Probabilistic Neural Network (2)
%    Bcl_qda               + Quadratic Discriminant Analysis 
%    Bcl_SMOTEBoost        + Smote Boost
%    Bcl_structure         + Classification using a Balu's structure
%    Bcl_svm               + Support Vector Machine  (4)
%    Bcl_svmplus           + SVM for 2 or more classes (4)
%    Bcl_tree              + Decision tree
%    Bcl_weakc             + Weak Classifier (3)
% 
% Feature analysis (fa)
%    Bfa_bestcorr          + Search best linear correlation 
%    Bfa_bestcorrn         + Search best linear combination correlation
%    Bfa_corrsearch        + Error estimation of linear or quadratic model
%    Bfa_dXi2              + Xi^2 distance between two vectors
%    Bfa_gmean             + Geometric mean
%    Bfa_jfisher           + score using Fisher discriminant
%    Bfa_kde               + Kernel density estimator 
%    Bfa_kde2d             + Kernel density estimator for 2D
%    Bfa_miparzen2         + Mutual Information using Parzen window
%    Bfa_score             + score of a set of features
%    Bfa_sp100             + score using sp=100%
%    Bfa_sqcorrcoef        + Squared-correlation coefficient
%    Bfa_vecsimilarity     + Normalized scalar product 
%
% Clustering (ct)
%    Bct_kmeans            + Clustering using k-means
%    Bct_knngraph2d        + KNN graph for 2D
%    Bct_meanshift         + Clustering using meanshift
%    Bct_medoidshift       + Clustering using medoidshif
%    Bct_neighbor          + Neighbor clustering
%    Bct_neighbor2D        + Neighbor clustering in 2D
%    Bct_spectralct        + Spectral clustering
%
% Performance evaluation (ev)
%    Bev_bootstrap         + Bootstrap evaluation
%    Bev_bootstrap0632     + Bootstrap 0.632 evaluation
%    Bev_confidence        + Confidence interval
%    Bev_confusion         + Confusion Matrix
%    Bev_crossval          + Cross-validation
%    Bev_holdout           + Holdout evaluation
%    Bev_jackknife         + Jackknife
%    Bev_performance       + Performance evaluation
%    Bev_roc               + ROC curve
%
% Multi-view analysis (mv)
%    Bmv_antisimetric      + Antisimetric matrix
%    Bmv_bundleafin        + Afin bundle adjustment
%    Bmv_bundleproj        + Projective bundle adjustment
%    Bmv_epidist           + Epipolar distance
%    Bmv_epiplot           + Plot of epipolar line
%    Bmv_epipoles          + Epipoles from fundamental matrix
%    Bmv_fundamental       + Fundamental matrix F from 2 projection matrices
%    Bmv_fundamentalRANSAC + RANSAC estimation of F from projection points
%    Bmv_fundamentalSIFT   + Estimation of F using SIFT points (1)
%    Bmv_fundamentalSVD    + SVD estimation of F from projection points
%    Bmv_guiproy2D         + GUI for projective transformation in 2D
%    Bmv_homographyRANSAC  + RANSAC estimation of homography matrix H
%    Bmv_homographySIFT    + Homography of two images using SIFT points (1)
%    Bmv_homographySVD     + SVD estimation of homography matrix H
%    Bmv_line2img          + Line to image conversion
%    Bmv_lines2point       + Intersction of two 2D lines
%    Bmv_matchSIFT         + Matching points between two images
%    Bmv_matrixp           + Perspective projection matrix 3D->2D
%    Bmv_matrixr2d         + Rotation matrix in 2D
%    Bmv_matrixr3d         + Rotation matrix in 3D
%    Bmv_points2line       + 2D line l that contains two 2D points 
%    Bmv_projective2D      + 2D projective transformation
%    Bmv_reco3d2           + 3D reconstruction in 2 views
%    Bmv_reco3dn           + 3D reconstruction in n views
%    Bmv_reco3dna          + 3D affine reconstruction in n views
%    Bmv_reproj3           + Reprojection of third view
%    Bmv_tqsift            + Target - query search using SIFT
%    Bmv_trifocal          + Trifocal tensors
%    Bmv_trifocalSVD       + Estimation of Trifocal Tensors using SVD decomposition
%
% Sequence processing (sq)
%    Bsq_des               + Description of a sequence (1)
%    Bsq_fundamental       + Fundamental matrices of a sequence (calib)
%    Bsq_fundamentalSIFT   + Fundamental matrix between two sequence views
%    Bsq_load              + Load an image sequence
%    Bsq_movie             + Movie of an image sequence
%    Bsq_multifundamental  + Fundamental matrices of a sequence (no calib)
%    Bsq_patch             + Patches of a sequence
%    Bsq_show              + Display of a sequence
%    Bsq_sort              + Sort of a sequence (1)
%    Bsq_stoplist          + Stoplist of a visual vocabulary
%    Bsq_stopout           + Filtering of stop words
%    Bsq_trifocal          + Trifocal tensors of a sequence
%    Bsq_vgoogle           + Search similar sequence images (1)
%    Bsq_visualvoc         + Visual vocabulary
%    Bsq_vocabulary        + Visual vocabulary of a sequence (1)
%
% Tracking (tr)
%    Btr_2                 + Matching in 2 views with epipolar constraint
%    Btr_3                 + Matching in 3 views with trifocal constraint
%    Btr_analysis          + Track analysis
%    Btr_classify          + Track classify using Multi-view windows
%    Btr_demo              + Track demo: detection using multi-view
%    Btr_detection         + Detection by tracking
%    Btr_gui               + Graphic user interface for tracking algorithm
%    Btr_join              + Join of matching points
%    Btr_merge             + Merge tracks with common matching points
%    Btr_plot              + Plot of tracks
%    Btr_reco3d            + 3D reconstruction of a track
%    Btr_sfm               + Structure from Motion
%    Btr_sfseq             + Structure from an image sequence
%    Btr_sift2             + Matching keypoints in 2 views of a sequence
%    Btr_siftn             + Matching keypoints in n views of a sequence
%    Btr_windows           + Multi-view windows (MVW)
%
% Miscellaneous
%    andsift               + SIFT descriptors of a region
%    clb                   + clt and delete(Bio_statusbar)
%    clt                   + close all and clear all
%    compare               + comparison between two variables
%    d2Y                   + conversion from label vector to binary matrix
%    distxy                + Euclidean distance between vectors of matrices
%    eigsort               + sorted eigen function
%    enterpause            + pause with "press enter to continue..."
%    i2h                   + inhomogeneous coordinates to homogeneous
%    howis                 + attributes of a matlab variable
%    h2i                   + homogeneous coordinates to inhomogeneous
%    imi                   + display image I
%    imshowc               + close all windows before display an image
%    imshows               + imshows(I) corresponds to imshow(I,[])
%    indices               + indices of a vector
%    max2                  + maximum of a matrix
%    mean2                 + average of a matrix
%    min2                  + minimum of a matrix
%    num2fixstr            + num2str using filling with '0' from left
%    posrandom             + interchange columns of a matrix randomly
%    sqdif                 + mean of Euclidean difference
%    sum2                  + sum off all elements of a matrix
%
% Notes:
% (1) It requires VLFeat Toolbox (see www.vlfeat.org).
% (2) It requires Neural Network Toolbox (see www.mathworks.com).
% (3) It requires Image Processing Toolbox (see www.mathworks.com).
% (4) It requires Bioinformatics Toolbox (see www.mathworks.com).
% (5) Implementation based on NetLab Toolbox (www.ncrg.aston.ac.uk/netlab).
% (6) Based on implementation by Heikkila & Ahonen 
%     (from http://www.cse.oulu.fi/MVG/Research/LBP).
% (7) Based on implementation by Anna Bosch 
%     (from http://www.robots.ox.ac.uk/~vgg/research/caltech/phog.html).
% (8) Based on implementation nipals.m by Geladi
%     (from http://www.cdpcenter.org/files/plsr).
% (9) It requires LIBSVM (see http://www.csie.ntu.edu.tw/~cjlin/libsvm/).
%(10) Copyright (c) 2004, Marcel Leutenegger.
