function Js = Bfa_score(X,d,options)

b      = options.b;
method = b.name;


switch lower(method)
    case 'mi' % mutual information
        if ~isfield(b,'param')
            dn = max(d)-min(d)+1; % number of classes

            p = ones(dn,1)/dn;
        else
            p = b.param;
        end
        Js = Bfa_mutualinfo(X,d,p);        
    case 'mr' % maximal relevance
        if ~isfield(b,'param')
            dn = max(d)-min(d)+1; % number of classes

            p = ones(dn,1)/dn;
        else
            p = b.param;
        end
        Js = Bfa_relevance(X,d,p);        
    case 'mrmr' % minimal redundancy and maximal relevance
        if ~isfield(b,'param')
            dn = max(d)-min(d)+1; % number of classes

            p = ones(dn,1)/dn;
        else
            p = b.param;
        end
        Js = Bfa_mRMR(X,d,p);        
    case 'fisher'
        if ~isfield(b,'param')
            dn = max(d)-min(d)+1; % number of classes

            p = ones(dn,1)/dn;
        else
            p = b.param;
        end
        Js = Bfa_jfisher(X,d,p);
    case 'sp100'
        Js = Bfa_sp100(X,d);
    otherwise
        ds = Bcl_structure(X,d,X,b);
        Js = Bev_performance(d,ds);
end
