% ds      = Bcl_ann(X,d,Xt,[])  Training & Testing together
% options = Bcl_ann(X,d,[])     Training only
% ds      = Bcl_ann(Xt,options) Testing only
%
% Toolbox: Balu
%    Simple Neural Network using Neural Network Toolbox of Matlab using
%    softmax
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.hidden is the number of hidden layers used in the IRLS algorithm
%       (default=10).
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.net contains information about the neural network
%       options.dmin contains min(d).
%       options.string is a 8 character string that describes the performed
%       classification (e.g., 'ann(10) ' means neural network with 10 hidden layers).
%
%    Example: Training & Test together:
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.hidden = 12;
%       ds = Bcl_ann(X,d,Xt,op);   % logistic - neural network
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.hidden = 12;
%       op = Bcl_ann(X,d,op);      % logistic - neural network - training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_ann(Xt,op);       % logistic - neural network - testing
%       p = Bev_performance(ds,dt) % performance on test data
%
% D.Mery, PUC-DCC, 2016
% http://dmery.ing.puc.cl

function [ds,options] = Bcl_ann(varargin)
[trainx,testx,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = sprintf('ann(%d) ',options.hidden);
if trainx
    dmin = min(d);
    dmax = max(d);
    d1   = d-dmin+1;
    k    = dmax-dmin+1;
    id   = eye(k);
    t    = (id(d1,:))';
    
    net = patternnet(options.hidden);
    net = train(net,X',t);
    
    options.net = net;
    options.dmin = dmin;
    ds = options;
    
    
end

if testx
    net = options.net;
    y = net(Xt');
    ds = (vec2ind(y)+options.dmin-1)';
end

end


