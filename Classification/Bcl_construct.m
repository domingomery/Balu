% This function is not a classifier!!!
% This function is called by Balu classifier functions (such as Bcl_lda) to
% build the training and testing data.

function [train,test,X,d,Xt,options] = Bcl_construct(varargin)

train = 0;
test  = 0;
switch nargin
    case 2 % [ds,options] = Bcl_svm(Xt,options)        % testing only
        X       = [];
        d       = [];
        Xt      = varargin{1};
        options = varargin{2};
        test    = 1;
    case 3   % options = Bcl_svm(X,d,options)          % training only
        X       = varargin{1};
        d       = double(varargin{2});
        if size(X,1)~=length(d)
            error('Length of label vector does not match number of instances.');
        end
        Xt      = [];
        options = varargin{3};
        train   = 1;
    case 4   % [ds,options] = Bcl_svm(X,d,Xt,options)  % training & test
        X       = varargin{1};
        d       = double(varargin{2});
        if size(X,1)~=length(d)
            error('Length of label vector does not match number of instances.');
        end
        Xt      = varargin{3};
        options = varargin{4};
        test    = 1;
        train   = 1;
    otherwise
        error('Bcl_construct: number of input arguments must be 2, 3 or 4.');
end
