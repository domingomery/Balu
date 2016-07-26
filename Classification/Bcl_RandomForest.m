% ds      = Bcl_RandomForest(X,d,Xt,[])  Training & Testing together
% options = Bcl_RandomForest(X,d,[])     Training only
% ds      = Bcl_RandomForest(Xt,options) Testing only
%
% Toolbox: Balu
%    Classifier using Random Forest. This implementation uses command
%    TreeBagger of Statistics and Machine Learning Toolbox of Matworks
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.NTrees number of trees
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'RForest ').
%
%    Example: Training & Test together:
%       load datagauss                     % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)              % plot feature space
%       op.NTrees  = 50;                   % Number of trees
%       ds = Bcl_RandomForest(X,d,Xt,op);  % RandomForest classifier
%       p = Bev_performance(ds,dt)         % performance on test data
%
%    Example: Training only
%       load datagauss                     % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)              % plot feature space
%       op.NTrees  = 50;                   % Number of trees
%       op = Bcl_RandomForest(X,d,op);     % RandomForest classifier - training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_RandomForest(Xt,op);      % RandomForest classifier - testing
%       p = Bev_performance(ds,dt)         % performance on test data
%
%    See also TreeBagger.
%
% D.Mery, PUC-DCC, 2015
% http://dmery.ing.puc.cl

function [ds,options] = Bcl_RandomForest(varargin)
[train,test,X,d,Xt,options] = Xconstruct(varargin{:});

options.string = 'RForest ';
if train
    options.B = TreeBagger(options.NTrees,X,d,'OOBPred','On','NVarToSample','all');
    ds = options;
end
if test
    dc = predict(options.B,Xt);
    ds = str2num(cell2mat(dc));
end

