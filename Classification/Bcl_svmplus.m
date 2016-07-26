% ds      = Bcl_svmplus(X,d,Xt,options)  Training & Testing together
% options = Bcl_svmplus(X,d,options)     Training only
% ds      = Bcl_svmplus(Xt,options)      Testing only
%
% Toolbox: Balu
%    Classifier using Support Vector Machine approach using Bioinformatics
%    Toolbox of Matlab using tree algorithm when number of classes is
%    2 or more. For two classes Bcl_svmplus calls Bcl_svm. For more than 2
%    Bsvmplus calls Btree.
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%
%       See help Bcl_svm to see options of Bcl_svmplus
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.string is a 8 character string that describes the performed
%       classification (e.g., 'svm+,4  ' means rbf-SVMplus).
%       options.svm = 1 means Bsvm was used.
%       options.svm = 0 means Btree with Bsvm was used.
%       See help Btree to see output options of Bsvmplus.
%
%    Example: Training & Test together:
%       [X,d] = Bds_gaussgen([2 1;1 2;2 2;1 1],ones(4,2)/4,500*ones(4,1));
%       [Xt,dt] = Bds_gaussgen([2 1;1 2;2 2;1 1],ones(4,2)/4,500*ones(4,1));
%       Bio_plotfeatures(X,d)        % plot feature space
%       op.kernel = 4;               % rbf
%       ds = Bcl_svmplus(X,d,Xt,op); % rbf-SVM
%       p = Bev_performance(ds,dt)   % performance on test data
%
%    Example: Training only
%       [X,d] = Bds_gaussgen([2 1;1 2;2 2;1 1],ones(4,2)/4,500*ones(4,1));
%       [Xt,dt] = Bds_gaussgen([2 1;1 2;2 2;1 1],ones(4,2)/4,500*ones(4,1));
%       Bio_plotfeatures(X,d)        % plot feature space
%       op.kernel = 4;               % rbf
%       op = Bcl_svmplus(X,d,op);    % rbf-SVM
%
%    Example: Testing only (after training only example):
%       ds = Bcl_svmplus(Xt,op);     % rbf-SVM classifier testing
%       p = Bev_performance(ds,dt)   % performance on test data
%
%    See also Bcl_svm, Bcl_tree, Bcl_libsvm.
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


function [ds,options] = Bcl_svmplus(varargin)

[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'svmplus ';

if ~isfield(options,'libsvm')
    options.libsvm = 0;
end


if train
    dmin = min(d);
    dmax = max(d);
    d    = d-dmin+1;

    if (dmax-dmin)==1
        op = options;
        if option.libsvm == 1
        options = Bcl_libsvm(X,d,op);
        options.svm = 1;
        else
        options = Bcl_svm(X,d,op);
        options.svm = 1;
        end
    else
        if options.libsvm == 1
            b.name = 'Bcl_libsvm';
        else
            b.name = 'Bcl_svm';
        end
        b.options.kernel = options.kernel;
        op = b;
        options = Bcl_tree(X,d,op);
        options.svm = 0;
        options.kernel = b.options.kernel;
    end    
    options.string = sprintf('svm+,%d  ',options.kernel);
    options.dmin = dmin;
    ds = options;
end
if test
    if options.svm
        if options.libsvm == 1
            ds = Bcl_libsvm(Xt,options) + options.dmin - 1;
        else
            ds = Bcl_svm(Xt,options) + options.dmin - 1;
        end
    else
        ds = Bcl_tree(Xt,options) + options.dmin - 1;
    end
end


