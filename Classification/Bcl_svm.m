% ds      = Bcl_svm(X,d,Xt,options)  Training & Testing together
% options = Bcl_svm(X,d,options)     Training only
% ds      = Bcl_svm(Xt,options)      Testing only
%
% Toolbox: Balu
%    Support Vector Machine approach using the Bioinformatics Toolbox.
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%
%       options.kernel defines the SVM-kernel as follows:
%
%       1: 'linear'      Linear kernel or dot product (default)
%       2: 'quadratic'   Quadratic kernel
%       3: 'polynomial'  Polynomial kernel (default order 3)
%       4: 'rbf'         Gaussian Radial Basis Function kernel
%       5: 'mlp'         Multilayer Perceptron kernel (default scale 1)
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
%       op.kernel = 4;
%       ds = Bcl_svm(X,d,Xt,op);   % rbf-SVM classifier
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.kernel = 4;
%       op = Bcl_svm(X,d,op);      % rbf-SVM classifier training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_svm(Xt,op);       % rbf-SVM classifier testing
%       p = Bev_performance(ds,dt) % performance on test data
%
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


function [ds,options] = Bcl_svm(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = sprintf('svm,%d   ',options.kernel);
if train
    k = {'linear','quadratic','polynomial','rbf','mlp'};
    options.svmStruct = svmtrain(X,d,'kernel_function',...
        cell2mat(k(options.kernel)),'method','ls');
    ds = options;
end
if test
    ds = svmclassify(options.svmStruct,Xt);
end


