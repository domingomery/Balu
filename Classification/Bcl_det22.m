% ds      = Bcl_det22(X,d,Xt,[])  Training & Testing together
% options = Bcl_det22(X,d,[])     Training only
% ds      = Bcl_det22(Xt,options) Testing only
%
% Toolbox: Balu
%    Quadratic Detector Design for ONLY two classes and two features
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.method is one of the following numbers:
%          1: Fisher
%          2: Best point of ROC curve
%          3: (1-Sp)->min @ Sn = 1
%          4: Sn->max @ Sp = 1
%          5: Neyman-Pearson: (1-Sp)-> min @ Sn = 90%
%          6: Least Squares Regression
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.dmin contains min(d).
%       options.dmax contains max(d).
%       options.w is the classifications vector.
%          w is the classification vector:
%             class 1: if x(1)*w(1)+...x(3) < 0
%             class 0: else
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'det21   ').
%
%    Example: Training & Test together:
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.method = 2;
%       ds = Bcl_det22(X,d,Xt,op); % Bdet22 with best ROC point
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.method = 2;
%       op = Bcl_det22(X,d,op);    % Bcl_det22 with best ROC point
%
%    Example: Testing only (after training only example):
%       ds = Bcl_det22(Xt,op);     % Bcl_det22 with best ROC point
%       p = Bev_performance(ds,dt) % performance on test data
%
%  See also Bdet21.
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [ds,options] = Bcl_det22(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});

options.string = 'det22   ';
if train
    if size(X,2)~=2
        error('Bcl_det22 works only for two features');
    end

    method = options.method;


    % Searching initial guest for w
    d0   = d;
    dmin = min(d0);
    dmax = max(d0);

    d(d0==dmin) = 0;
    d(d0==dmax) = 1;

    m1 = mean(X(d0==dmax,:));
    x0 = m1(1);
    y0 = m1(2);
    r = 0.3;
    wini = [1 0 1 -2*x0 -2*y0 (x0^2+y0^2-r^2)]'; %[a b c d e f], pp 9 Hartley

    if method == 1
        wini = fminsearch(@Bdeterr22,wini,[],X,d,2);
    end
    w = fminsearch(@Bdeterr22,wini,[],X,d,method);
    options.w = w;
    options.dmin = dmin;
    options.dmax = dmax;
    ds = options;
end
if test
    XXt = [Xt(:,1).*Xt(:,1) Xt(:,1).*Xt(:,2) Xt(:,2).*Xt(:,2) Xt(:,1) Xt(:,2) ones(size(Xt,1),1)];
    ds0 = double(XXt*options.w < 0); % Linear Classification
    ds  = zeros(size(ds0));
    ds(ds0==0) = options.dmin;
    ds(ds0==1) = options.dmax;
end
end


function [e,C]  = Bdeterr22(w,X,d,method)

XX = [X(:,1).*X(:,1) X(:,1).*X(:,2) X(:,2).*X(:,2) X(:,1) X(:,2) ones(size(X,1),1)];
ds = XX*w < 0; % linear classification

TP = sum(ds.*d);
TN = sum(not(ds).*not(d));
FN = sum(not(ds).*d);
FP = sum(ds.*not(d));

Sn  = TP/(TP+FN); % Sensitibity
Sp1 = FP/(TN+FP); % 1 - Specificity

lambda = 1e6;

switch method
    case 1 % Fisher, i.e., trace(inv(Cw)*Cb)->max
        z0 = X(ds==0,:);
        z1 = X(ds==1,:);
        zb0 = mean(z0)';
        zb1 = mean(z1)';
        Cw = cov(z0) + cov(z1);
        Cb = zb0*zb0'+zb1*zb1';
        J = trace(inv(Cw)*(Cb)); %Fisher discriminant
        e = -J;
    case 2 % Best point of ROC, i.e., (1-Sp)^2+(1-Sn)^2->Min
        e = sqrt(Sp1^2+(1-Sn)^2);
    case 3 % FP->min @ Sn = 100%
        % e = lambda*FN+FP;
        e = lambda*abs(Sn-1)+FP;
    case 4 % FN->min @ Sp = 100%
        e = FN+lambda*FP;
    case 5 % Neyman Pearson: FP->min @Sn = 90%
        e = lambda*abs(Sn-0.9)+FP;
    case 6 % Least Squares: norm(Xw-d)->min
        d5 = -2*(d-0.5);
        e = norm(d5-XX*w);
    otherwise
        error('Bdet22: method not defined.');

end

C = [TP TN FP FN];

% disp(sprintf('TP=%3d TN=%3d FP=%3d FN=%3d e=%7.4f',TP,TN,FP,FN,e))

end