% selec = Bfsransac(X,d,m,show,method,param,param2,param3)
%
% Toolbox: Balu
%    Sequential Forward Selection for fatures X according to ideal
%    classification d. m features will be selected.
%    method = 'fisher' uses Fisher objetctive function (in this case param
%    is the a priori probability of each class, if not given it will be
%    assumed constant).
%    method = 'sp100' uses as criteria Sp @Sn=100%.
%    method could be any classifier implemented in Balu with the
%    corresponding parameters, e.g., method = 'knn' and param = 10, it will
%    select the best m features for 10 nearest neighbours.
%    show = 1 display results.
%    selec is the indices of the selected features.
%
% D.Mery, PUC-DCC, Ago. 2010
% http://dmery.ing.puc.cl
%

function selec = Bfsransac(X,d,m,show,method,varargin)

s0 = Bfsclean(X);


f = X(:,s0);
nf = size(f,2);

selec  = []; %selected features
Jmax = 0;

dn = max(d)-min(d)+1; % number of classes

if not(exist('show','var'))
    show = 0;
end

if not(exist('method','var'))
    method = 'fisher';
end


k = 0;
while (k<=300)
    [f1,d1] = Bstratify(f,d,0.67);
    s = Bfsfs(f1,d1,m);
    fs = f(:,s);
    
    
    switch lower(method)
        case 'fisher'
            if (not(exist('param','var')))
                p = ones(dn,1)/dn;
            end
            Js = jfisher(fs,d,p);
        case 'sp100'
            Js = sp100(fs,d);
        otherwise
            e = ['ds = ' method '(fs,d,fs,varargin{:});'];
            eval(e);
            Js = Bperformance(d,ds);
    end
    if (Js>Jmax)
        selec = s;
        Jmax = Js;
        if show
            fprintf('Jmax = %8.4f\n',Jmax);
        end
    end
    k = k+1;
end
