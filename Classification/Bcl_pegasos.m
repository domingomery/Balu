% ds      = Bcl_pegasos(X,d,Xt,[])  Training & Testing together
% options = Bcl_pegasos(X,d,[])     Training only
% ds      = Bcl_pegasos(Xt,options) Testing only
%
% Toolbox: Balu
%    Classifier using Pegasos Support Vector Machine approach using
%    VLFeat Toolbox of Matlab.
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.lambda is the regularization parameter, can be scalar or a vector
%          if lambda is a vector, Bpegasos test for every elemente of
%          lambda the best performance in the reclassification.
%       options.n is the number of times it will be tested for each lambda.
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.dmin contains min(d).
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'pegasos ').
%
%    Example using linear SVM
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.lambda = 0.1:0.1:1.2;
%       op.n = 10;
%       ds = Bcl_pegasos(X,d,Xt,op);
%       p = Bev_performance(dt,ds)
%
%    Example using X^2 SVM
%       PX  = vl_homkermap(X',1,.6,'kchi2')';
%       PXt = vl_homkermap(Xt',1,.6,'kchi2')';
%       ds = Bcl_pegasos(PX,d,PXt,op);
%       p = Bev_performance(dt,ds)
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.lambda = 0.1:0.1:1.2;
%       op.n = 10;
%       op = Bcl_pegasos(X,d,op);  % pegasos - training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_pegasos(Xt,op);   % pegasos - testing
%       p = Bev_performance(ds,dt) % performance on test data
%
% D.Mery, PUC-DCC, Jun. 2010
% http://dmery.ing.puc.cl




function ds = Bcl_pegasos(varargin)

if ~exist('vl_pegasos','file')
    error('Bpegasos: This function requires the VLFeat Toolbox.');
end

[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'pegasos ';
if train
    lambda = options.lambda;
    n      = options.n;

    dmin = min(d);
    dmax = max(d);
    c = dmax-dmin+1;
    if c~=2
        error('Bcl_pegasos is implemented for only two classes.')
    end
    d=d-dmin;
    pmax = 0;
    for l = lambda
        for i=1:n
            w = vl_pegasos(X',int8(d),l);
            dt=double(X*w>0.5);
            p = Bev_performance(d,dt);
            if p>pmax
                pmax = p;
                % lmax = lambda;
                options.wmax = w;
            end
        end
    end
    options.dmin = dmin;
    ds = options;
end
if test
    p  = Xt*options.wmax-0.5;
    sc = ones(size(p))./(abs(p)+1e-5);
    ds = double(p>0)+options.dmin;
    ds = Bcl_outscore(ds,sc,options);
end

