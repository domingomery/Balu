% ds      = Bcl_det21(X,d,Xt,[])  Training & Testing together
% options = Bcl_det21(X,d,[])     Training only
% ds      = Bcl_det21(Xt,options) Testing only
%
% Toolbox: Balu
%    Linear Detector Design for ONLY two classes and two features
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.method is one of the following numbers:
%          0: Fisher0
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
%       ds = Bcl_det21(X,d,Xt,op); % Bdet21 with best ROC point
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.method = 2;
%       op = Bcl_det21(X,d,op);    % Bdet21 with best ROC point
%
%    Example: Testing only (after training only example):
%       ds = Bcl_det21(Xt,op);     % Bdet21 with best ROC point
%       p = Bev_performance(ds,dt) % performance on test data
%
% See also Bcl_det22.
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


function [ds,options] = Bcl_det21(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});

options.string = 'det21   ';
if train
    if size(X,2)~=2
        error('Bdet21 works only for two features');
    end

    method = options.method;


    % Searching initial guest for w
    d0   = d;
    dmin = min(d0);
    dmax = max(d0);

    i0 = find(d0==dmin);
    i1 = find(d0==dmax);
    d(i0) = 0;
    d(i1) = 1;
    n0 = length(i0);
    n1 = length(i1);
    n = n0+n1;

    % means
    m0 = mean(X(i0,:))';
    m1 = mean(X(i1,:))';

    %covariances
    S0 = cov(X(i0,:));
    S1 = cov(X(i1,:));

    %within covariance
    Sw = (n0*S0+n1*S1)/(n-2);

    w = inv(Sw)*(m0-m1);
    p1 = 0.5;
    p0 = 0.5;
    w0 = -0.5*(m0+m1)'*inv(Sw)*(m0-m1)-log(p1/p0);
    wfis = [w;w0];

    % least squares
    XX = [X ones(size(X,1),1)];
    d2 = max(d);
    d1 = min(d);
    d = (d-d1)*1/(d2-d1);
    d5 = -2*(d-0.5);
    wls = inv(XX'*XX)*XX'*d5;

    switch method
        case 0
            w = wfis;
        case {1,2,3,4,5}
            if method>1
                emin = Inf;
                ax = axis;
                x0 = ax(1);
                dx = ax(2)-ax(1);
                y0 = ax(3);
                dy = ax(4)-ax(3);

                for kk=1:10000
                    nn = [x0 + rand*dx;y0 + rand*dy; 1];
                    mm = [x0 + rand*dx;y0 + rand*dy; 1];
                    w = cross(nn,mm);
                    e = Bdeterr21(w,X,d,method);
                    if e<emin
                        emin = e;
                        wini = w;
                    end
                end
            else
                wini = wfis;
            end
            w = fminsearch(@Bdeterr21,wini,[],X,d,method);
        case 6
            w = wls;
        otherwise
            error('Bdet21: method not defined.');
    end
    options.w = w;
    options.dmin = dmin;
    options.dmax = dmax;
    ds = options;
end
if test
    XXt = [Xt ones(size(Xt,1),1)];
    ds0 = double((XXt*options.w) < 0); % Linear Classification
    ds  = zeros(size(ds0));
    ds(ds0==0) = options.dmin;
    ds(ds0==1) = options.dmax;
end
end


function [e,C]  = Bdeterr21(w,X,d,method)

XX = [X ones(size(X,1),1)];
ds = XX*w < 0; % Linear Classification

TP = sum(ds.*d);
TN = sum(not(ds).*not(d));
FN = sum(not(ds).*d);
FP = sum(ds.*not(d));

Sn  = TP/(TP+FN); % Sensitibity
Sp1 = FP/(TN+FP); % 1 - Specificity

lambda = 1e5;

switch method
    case 0
        e = NaN; %not defined
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
end

C = [TP TN FP FN];

end
