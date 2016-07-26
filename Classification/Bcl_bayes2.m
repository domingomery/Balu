% ds      = Bcl_bayes2(X,d,Xt,options)  Training & Testing together
% options = Bcl_bayes2(X,d,options)     Training only
% ds      = Bcl_bayes2(Xt,options)      Testing only
%
% Toolbox: Balu
%    Bayes classifier for ONLY two features and two classes
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.p is the prior probability, if p is not given, it will be estimated
%         proportional to the number of samples of each class.
%       options.show = 1 displays results.
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.dmin contains min(d).
%       options.string is a 8 character string that describes the performed
%       classification (in this case 'bayes2  ').
%       options.m1,m2,b1,b2,C are parameters of the classifier.
%
%
%    Example: Training & Test together:
%       load datagauss              % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)       % plot feature space
%       op.p = [0.25 0.75];         % prior probability
%       op.show = 1;                % display results
%       ds = Bcl_bayes2(X,d,Xt,op); % Bayes classifier
%       p = Bev_performance(ds,dt)  % performance on test data
%
%    Example: Training only
%       load datagauss              % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)       % plot feature space
%       op.p = [0.25 0.75];         % prior probability
%       op.show = 1;                % display results
%       op = Bcl_bayes2(X,d,op);    % Bayes classifier
%
%    Example: Testing only (after training only example):
%       ds = Bcl_bayes2(Xt,op);     % Bayes classifier
%       p = Bev_performance(ds,dt)  % performance on test data
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


function [ds,options] = Bcl_bayes2(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'bayes2  ';

if train

    dmin = min(d);
    
    d = d - dmin + 1;


    p = options.p;
    show = options.show;

    col = 'gbrcmykb';

    n = size(X,1);

    N = max(d);
    if show
        figure(1)
        Bio_plotfeatures(X,d)
    end
    x0 = min(X(:,1));
    x1 = max(X(:,1));

    y0 = min(X(:,2));
    y1 = max(X(:,2));

    NN = 100;
    dx = (x1-x0)/NN;
    dy = (y1-y0)/NN;
    H = zeros(NN,NN,N);
    i0 = 0.500001;
    i1 = NN+1-i0;
    m1 = (i1-i0)/(x1-x0);
    b1 = i0 - m1*x0;
    m2 = (i1-i0)/(y1-y0);
    b2 = i0 - m2*y0;

    for k = 1:n
        i = fix(m1*X(k,1)+b1+0.5);
        j = fix(m2*X(k,2)+b2+0.5);
        H(i,j,d(k)) = H(i,j,d(k)) + 1;
    end
    G = zeros(size(H));
    for k=1:N
        H(:,:,k) = H(:,:,k)/sum(sum(H(:,:,k)))*p(k);
        G(:,:,k) = Bgaussestimator(H(:,:,k));
        if show
            figure(2)
            mesh(H(:,NN:-1:1,k))
            title('Histogram');
            hold on
            figure(3)
            mesh(G(:,NN:-1:1,k))
            title('Modeled gaussian pdf');
            hold on
            drawnow
        end
    end

    C = zeros(NN,NN);
    for i=1:NN
        for j=1:NN
            [c,r] = max(G(i,j,:));
            if (c>0)
                C(i,j) = r;
            end
        end
    end
    options.m1 = m1;
    options.m2 = m2;
    options.b1 = b1;
    options.b2 = b2;
    options.C = C;
    options.dmin = dmin;
    ds = options;

    if show
        figure(1)
        hold on
        for x = x0:dx:x1
            i = fix(m1*x+b1+0.5);
            for y = y0:dy:y1
                j = fix(m2*y+b2+0.5);
                plot(x,y,[col(C(i,j)+1) ':'])
            end
            % drawnow
        end
        title('Test Data')
    end
end
if test
    nt   = size(Xt,1);
    m1   = options.m1;
    m2   = options.m2;
    b1   = options.b1;
    b2   = options.b2;
    C    = options.C;
    dmin = options.dmin;
    NN = size(C,1);

    ds = zeros(nt,1);
    for k = 1:nt
        i = fix(m1*Xt(k,1)+b1+0.5);
        if (i>NN)
            i = NN;
        end
        if (i<1)
            i = 1;
        end
        j = fix(m2*Xt(k,2)+b2+0.5);
        if (j>NN)
            j = NN;
        end
        if (j<1)
            j = 1;
        end
        ds(k) = C(i,j)+dmin-1;
    end
end
end


function G = Bgaussestimator(H)
H = conv2(H,ones(5,5)/25,'same');
n = size(H);
if (length(n) == 2)
    if or(n(1)==1,n(2)==1) % one dimension

    else % two dimension
        [ii,jj] = find(H>0);
        im = mean(ii);
        jm = mean(jj);
        si = (im-min(ii))/2;
        sj = (jm-min(jj))/2;
        hmax = max(H(:));
        th0 = [hmax im jm si sj 0];
        th = fminsearch(@Bgausserr2,th0,[],H);
        [err,G]  = Bgausserr2(th,H);


    end
else
    error('Bbayes2 works with only two dimension');
end
end




function [err,G] = Bgausserr2(th,H)
% th = [hmax im jm si sj al];
hmax = th(1); % maximum value of Gauss function
im   = th(2); % mean value in i direction
jm   = th(3); % mean value in j direction
si   = th(4); % std in i direction
sj   = th(5); % std in j direction
al   = th(6); % orientation
n = size(H,1);

a = (cos(al)/sj)^2 + (sin(al)/si)^2;
b = -sin(2*al)/sj^2 + sin(2*al)/si^2 ;
c = (sin(al)/sj)^2 + (cos(al)/si)^2;

[j, i] = meshgrid(1:n,1:n);

G = hmax*exp( - (a*(j-jm).^2 + b*(j-jm).*(i-im) + c*(i-im).^2)) ;

D = abs(G(:)-H(:));
err = mean(D);
end


