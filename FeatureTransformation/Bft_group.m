function Y = Bft_group(X,method,K)

switch lower(method)
    case 'hist' 
        if ~exist('K','var')
            K = max(X);
        end
        Y            = hist(X,1:K);
    case 'histnorm' 
        if ~exist('K','var')
            K = max(X);
        end
        Y            = hist(X,1:K);
        Y            = Y/sum(Y);
    case 'homkermap'
        if ~exist('K','var')
            K = max(X);
        end
        h    = hist(X,1:K);
        h    = h/sum(h);
        Y    = vl_homkermap(h', 1, 'kchi2', 'gamma', .5) ;
    case 'sum'
        Y    = sum(X);
    otherwise
        error('%d not defined.\n',method)

end
