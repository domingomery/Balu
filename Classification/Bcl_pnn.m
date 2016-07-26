% ds      = Bcl_pnn(X,d,Xt,options)  Training & Testing together
% options = Bcl_pnn(X,d,options)     Training only
% ds      = Bcl_pnn(Xt,options)      Testing only
%
% Toolbox: Balu
%    Probabilistic neural network (Neural Network Toolbox required).
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options must be [] (it is for future purposes)
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.dmin contains min(d).
%       options.net contains information about the neural network
%       (from function newpnn from Nueral Network Toolbox).
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'pnn     ').
%
%    Example: Training & Test together:
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       ds = Bcl_pnn(X,d,Xt,[]);   % pnn classifier
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op = Bcl_pnn(X,d,[]);      % pnn - training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_pnn(Xt,op);       % pnn - testing
%       p = Bev_performance(ds,dt) % performance on test data
%
% (c) D.Mery, PUC-DCC, May 2010
% http://dmery.ing.puc.cl


function [ds,options] = Bcl_pnn(varargin)
if isempty(ver('nnet'))
    error('Bcl_pnn requires the Neural Network Toolbox.');
end
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'pnn     ';
if train
    dmin = min(d);
    options.dmin = dmin;
    d    = d-dmin+1;
    T    = ind2vec(d');
    options.net  = newpnn(X',T,1);
    ds = options;
end
if test
    A    = sim(options.net,Xt');
    ds   = vec2ind(A)'+options.dmin-1;
end
