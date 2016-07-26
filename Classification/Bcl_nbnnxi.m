% ds =  Bcl_nbnnxi(X,d,Xt,D);
%
% Toolbox: Balu
%    Naive Bayes Nearest Neighbor for histograms using Xi distance
%
%    Design data:
%    X is the feature matrix having M histograms of D bins each per sample
%    d is the ideal classification for X
%
%    Test data:
%    Xt is the feature matrix
%
%    options.D: number of bins of the histograms
%
%    ds is the classification of Xt according the designed classifier
%
% D.Mery, PUC-DCC, Mar 2010
% http://dmery.ing.puc.cl
%
function ds = Bcl_nbnnxi(varargin)

[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = 'nbnnxi  ';
if train
    options.X = X;
    options.d = d;
    ds = options;
end
if test
    D    = options.D;
    show = options.show;
    d    = options.d;
    X    = options.X;
    K    = max(d)-min(d)+1; % number of classes
    M    = size(X,2)/D;     % number of historgrams
    N    = size(X,1);
    Nt   = size(Xt,1);
    ds   = zeros(Nt,1);
    dmin = min(d);
    d    = d-dmin+1;
    if show
        ff = Bio_statusbar('Bcl_nbnnxi');
    end
    for x=1:Nt
        dd = zeros(M,K);
        for q=1:M
            q0 = (q-1)*D+1;
            q1 = q*D;
            XX = Xt(x,q0:q1);
            Xi = zeros(N,1);
            for y=1:N
                YY    = X(y,q0:q1);
                s     = XX+YY;
                d2    = (XX-YY).*(XX-YY);
                Xi(y) = sum(d2(s>0)./s(s>0));
            end
            for j=1:K
                dd(q,j) = min(Xi(d==j));
            end
        end
        dsum = sum(dd)';
        [i,j] = min(dsum);
        ds(x) = j;
        if show
            ff = Bio_statusbar(x/Nt,ff);
            fprintf('%d/%d -> class = %d \n',x,Nt,j);
        end
    end
    ds = ds+dmin-1;
    if show
        delete(ff);
    end
end
