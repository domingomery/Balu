% ds      = Bcl_svm(X,d,Xt,options)  Training & Testing together
% options = Bcl_svm(X,d,options)     Training only
% ds      = Bcl_svm(Xt,options)      Testing only
%
% Toolbox: Balu
%    Support Vector Machine approach using the LIBSVM(*).
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%
%       options.kernel is a string that defines the options of LIBSVM as 
%       follows:
%
%          -s svm_type : set type of SVM (default 0)
%                  0 -- C-SVC
%                  1 -- nu-SVC
%                  2 -- one-class SVM
%                  3 -- epsilon-SVR
%                  4 -- nu-SVR
%          -t kernel_type : set type of kernel function (default 2)
%                  0 -- linear: u'*v
%                  1 -- polynomial: (gamma*u'*v + coef0)^degree
%                  2 -- radial basis function: exp(-gamma*|u-v|^2)
%                  3 -- sigmoid: tanh(gamma*u'*v + coef0)
%          -d degree : set degree in kernel function (default 3)
%          -g gamma : set gamma in kernel function (default 1/num_features)
%          -r coef0 : set coef0 in kernel function (default 0)
%          -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
%          -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
%          -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
%          -m cachesize : set cache memory size in MB (default 100)
%          -e epsilon : set tolerance of termination criterion (default 0.001)
%          -h shrinking: whether to use the shrinking heuristics, 0 or 1 (default 1)
%          -b probability_estimates: whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
%          -wi weight: set the parameter C of class i to weight*C, for C-SVC (default 1)%
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.svmStruct contains information about the trained classifier
%       (from function svmtrain of Bioinformatics Toolbox).
%       options.string is a 8 character string that describes the performed
%       classification (e.g., 'Bsvm,4  ' means rbf-SVM).
%
%    Example: Training & Test together:
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.kernel = '-t 2';
%       ds = Bcl_libsvm(X,d,Xt,op);% rbf-SVM classifier
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.kernel = '-t 2';
%       op = Bcl_libsvm(X,d,op);   % rbf-SVM classifier training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_libsvm(Xt,op);    % rbf-SVM classifier testing
%       p = Bev_performance(ds,dt) % performance on test data
%
%  See also Bcl_svm, Bcl_svmplus.
% 
% Reference:
%    Chang, Chih-Chung and Lin, Chih-Jen (2011): LIBSVM: A library for 
%    support vector machines, ACM Transactions on Intelligent Systems and 
%    Technology, 2(3):27:1--27:27.
%    http://www.csie.ntu.edu.tw/~cjlin/libsvm
%
%    WARNING: When installing libsvm you must rename the compiled files
%    svmtrain and svmpredict (in libsvmxxx/matlab folder) as
%    libsvmtrain and libsvmpredict. Additionally, libsvmxxx/matlab folder  
%    must be included in the matlab path.
%
% D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [ds,options] = Bcl_libsvm(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string      = ['libsvm,%d' options.kernel];
if train
    s = options.kernel;
    options.model   = libsvmtrain(d,X,s);
    
    outp = 0;
    if isfield(options,'output')
       outp = options.output;
    end
    if outp == 3
       [ds,a,sc]     = libsvmpredict(ones(size(X,1),1),X,options.model);
       [ysc,param]   = Bft_sigmoid(sc,d);
       options.param = param;
    end
    ds              = options;
end
if test
    [ds,a,sc] = libsvmpredict(ones(size(Xt,1),1),Xt,options.model);
    % sc = abs(sc-ds);
    % sc = -sc;
    ds = Bcl_outscore(ds,sc,options);
end