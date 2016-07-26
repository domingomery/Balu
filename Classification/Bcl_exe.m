% ds      = Bcl_exe(bname,X,d,Xt,options)  Training & Testing together
% options = Bcl_exe(bname,X,d,options)     Training only
% ds      = Bcl_exe(bname,Xt,options)      Testing only
%
% Toolbox: Balu
%    Classification using Balu classifier bname and options.
%
%    Design data:
%       bname can be any name of a Balu classifier with or without 'B' (e.g. 
%          'Bsvm', 'Bknn', or 'svm', 'knn', etc.)
%       etc.)
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options are the options of the Balu classifier bname. (e.g., for
%       Bsvm see help Bsvm to see the corresponding options).
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data (one column per classifier)
%       options are the corresponding output options of the classifier.
%
%    Example: Training & Test together:
%       load datagauss                   % simulated data (2 classes, 2 features)
%       op.kernel = 4;
%       ds = Bcl_exe('svm',X,d,Xt,op);   % Bcl_svm classifier ('svm' or 'Bcl_svm')
%       p = Bev_performance(ds,dt)       % performance
%
%
%    Example: Training only
%       load datagauss                   % simulated data (2 classes, 2 features)
%       op.kernel = 4;
%       op = Bcl_exe('Bcl_svm',X,d,op);  % Bcl_svm classifier ('svm' or 'Bcl_svm')
%
%    Example: Testing only (after training only example):
%       ds = Bcl_exe('Bcl_svm',Xt,op);   % Bcl_svm classifier ('svm' or 'Bcl_svm')
%       p  = Bev_performance(ds,dt)
%
%    See also Bcl_structure.
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [ds,options] = Bcl_exe(bname,varargin)

if bname(1) == 'B'
   Bn = bname;
else
    Bn = ['Bcl_' bname];
end

str = [Bn(5:end) '        '];
options.string = str(1:8);

[ds,options] = feval(Bn,varargin{:});
