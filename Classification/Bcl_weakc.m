% ds      = Bcl_weakc(X,d,Xt,options)  Training & Testing together
% options = Bcl_weakc(X,d,options)     Training only
% ds      = Bcl_weakc(Xt,options)      Testing only
%
% Toolbox: Balu
%    Weak classifier for one feature X using Otsu method.
%    Design data:
%       X is a column vector with only one feature
%       d is the ideal classification for X
%       options must be [] (it is for future purposes)
%
%    Test data:
%       Xt is a column vector with only one feature
%
%    Output:
%       ds is the classification on test data
%       options.dmin contains min(d).
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'weakc   ').
%       options.th and options.p as follows:
%          p : 1 means ds = Xt<th
%             -1 means ds = Xt>th
%          th: threshold
%
%    Example: Training & Test together:
%       load datagauss                    % simulated data (2 classes, 2 features)
%       x  = X(:,1);                      % first feature
%       xt = Xt(:,1);                     % test data
%       Bio_plotfeatures(x,d)             % plot feature space
%       ds = Bcl_weakc(x,d,xt,[]);        % weak classifier
%       p = Bev_performance(ds,dt)        % performance on test data
%
%    Example: Training only
%       load datagauss                    % simulated data (2 classes, 2 features)
%       x  = X(:,1);                      % first feature
%       Bio_plotfeatures(x,d)             % plot feature space
%       op = Bcl_weakc(x,d,[]);           % weak classifier
%       hold on
%       plot([op.th op.th],[0 0.15],'k:') % threshold line
%
%    Example: Testing only (after training only example):
%       xt = Xt(:,1);                     % test data
%       ds = Bcl_weakc(xt,op);            % weak classifier
%       hold on
%       plot([op.th op.th],[0 0.15],'k:') % threshold line
%       p = Bev_performance(ds,dt)        % performance on test data
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [ds,options] = Bcl_weakc(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'weakc   ';
if train
    dmin = min(d);
    dmax = max(d);
    if (dmax-dmin~=1)
        error('Bcl_weakc works with only two classes.')
    end
    if size(X,2) > 1
        error('Bcl_weakc works with only one feature.')
    end    
    y = d-dmin;
    [xn,a,b]  = Bft_norm(X,0);
    th  = (graythresh(xn)-b)/a;
    ys1 = X<th;
    ys2 = not(ys1);
    e1  = sum(abs(ys1-y));
    e2  = sum(abs(ys2-y));
    if e1<e2
        p  = 1;
    else
        p  = -1;
    end
    options.dmin = dmin;
    options.th   = th;
    options.p    = p;
    ds           = options;
end
if test
    if size(Xt,2) > 1
        error('Bcl_weakc works with only one feature.')
    end    
    ys = options.p*Xt<(options.p*options.th);
    ds = ys+options.dmin;
end

