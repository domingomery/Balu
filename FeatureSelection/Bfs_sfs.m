% selec = Bfs_sfs(X,d,options)
%
% Toolbox: Balu
%    Sequential Forward Selection for fatures X according to ideal
%    classification d. optins.m features will be selected.
%    options.b.method = 'fisher' uses Fisher objetctive function.
%    options.b.method = 'sp100' uses as criteria Sp @Sn=100%.
%    options.b can be any Balu classifier structure (see example).
%    options.show = 1 display results.
%    selec is the indices of the selected features.
%
% Example 1: SFS using Fisher dicriminant
%    load datareal
%    op.m = 10;                     % 10 features will be selected
%    op.show = 1;                   % display results
%    op.b.name = 'fisher';          % SFS with Fisher
%    s = Bfs_sfs(f,d,op);           % index of selected features
%    X = f(:,s);                    % selected features
%    Xn = fn(s,:)                   % list of feature names
%    op_lda.p = [];
%    ds = Bcl_lda(X,d,X,op_lda);    % LDA classifier
%    p = Bev_performance(d,ds)      % performance with sfs 
%
% Example 2: SFS using a KNN classifier
%    load datareal
%    op.m = 10;                     % 10 features will be selected
%    op.show = 1;                   % display results
%    op.b.name = 'knn';             % SFS with KNN
%    op.b.options.k = 5;            % 5 neighbors
%    s = Bfs_sfs(f,d,op);           % index of selected features
%    X = f(:,s);                    % selected features
%    Xn = fn(s,:)                   % list of feature names
%
% Example 3: SFS using sp100 criterion
%    load datareal
%    op.m = 10;                     % 10 features will be selected
%    op.show = 1;                   % display results
%    op.b.name = 'sp100';           % SFS with sp100 criterion
%    s = Bfs_sfs(f,d,op);           % index of selected features
%    X = f(:,s);                    % selected features
%    Xn = fn(s,:)                   % list of feature names
%
% (c) D.Mery, PUC-DCC, Jul. 2011
% http://dmery.ing.puc.cl

function selec = Bfs_sfs(X,d,options)

m      = options.m;
show   = options.show;
if ~isfield(options,'force')
   options.force = 0;
end
force = options.force;


f      = X;
N      = size(f,2);
selec  = []; %selected features
J      = zeros(m,1);
% Jmax   = 0;
k      = 0;

if show
    ff = Bio_statusbar('SFS progress');
end
while (k<m)
    if show
        ff = Bio_statusbar(k/m,ff);
    end
    fnew = 0;
    Jmax = -Inf;
    for i=1:N
        if (k==0) || (sum(selec==i)==0)
            s = [selec; i];
            fs = f(:,s);
            Js = Bfa_score(fs,d,options);
            if (Js>=Jmax)
                ks = i;
                Jmax = Js;
                fnew = 1;
            end
        end
    end
    if (fnew)
        selec = [selec; ks];
        k = k + 1;
        J(k) = Jmax;
        if show
            clf
            bar(J)
            
            fprintf('Jmax = %8.4f\n',Jmax);
            hold on
            for i=1:length(selec)
                text(i-0.4,J(i)*1.05,sprintf('%d',selec(i)));
            end
            pause(0)
        end
    else
        disp('Bfs_sfs: no more improvement. Sequential search for feature selection is interrupted.');
        if and(force,(k<m))
           fprintf('Bfs_sfs: Warning! %d random features were selected in order to have\n',m-k);
           fprintf('                  %d selected features (options.force is true).\n',m);
           t = 1:N;
           t(selec) = [];
           n = length(t);
           x = rand(n,1);
           [i,j] = sort(x);
           selec = [selec; t(j(1:m-k))' ];
         end

        k = 1e10;
    end
end
if show
    delete(ff);
end
