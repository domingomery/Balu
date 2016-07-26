% [bcs,selec,sp] = Bcl_balu(X,d,bcl,bfs,options)
%
% Toolbox: Balu
%    Feature and classifier selection tool.
%    Exhaustive search of the best classifier of the classifiers given in
%    bcl structure using the features selected by feature selection
%    algorithms given in bfs structure.
%
%    X features
%    d ideal      classification
%    bcl          structure of classifiers
%    bfs          structure of feature selectors
%    options.Xn   name of the features (optional)
%    options.v    cross validation folder (default = 10)
%    options.c    cross validation confidence interval (default = 0.95)
%    options.m    maximum number of features to be selected
%    options.show display results (deafult = 1)
%
%    bcs: selected classifier
%    selec: selected features
%    sp is a Nc x m x Nq matrix with information about the performance
%    sp(i,j,q) contains cross validation information obtained by classifier
%    bcl(i), using the first j select features given by feature selection
%    algorithm q.
%
% Example:
% load datareal % f & fn is loaded
%
% bcl(1).name = 'knn';   bcl(1).options.k = 5;      % KNN with 5 neighbors
% bcl(2).name = 'knn';   bcl(2).options.k = 7;      % KNN with 7 neighbors
% bcl(3).name = 'knn';   bcl(3).options.k = 9;      % KNN with 9 neighbors
% bcl(4).name = 'lda';   bcl(4).options.p = [];     % LDA
% bcl(5).name = 'qda';   bcl(5).options.p = [];     % QDA
% bcl(6).name = 'svm';   bcl(6).options.kernel = 4; % rbf-SVM
% bcl(7).name = 'maha';  bcl(7).options = [];       % Euclidean distance
% bcl(8).name = 'dmin';  bcl(8).options = [];       % Mahalanobis distance
%
%
% bfs(1).name = 'sfs';   bfs(1).options.b.name    = 'fisher';
% bfs(2).name = 'sfs';   bfs(2).options.b.name    = 'knn'; bfs(2).options.b.options.k = 5;
% bfs(3).name = 'sfs';   bfs(3).options.b.name    = 'sp100';
% bfs(4).name = 'sfs';   bfs(4).options.b.name    = 'knn'; bfs(4).options.b.options.k = 7;
% bfs(5).name = 'rank';  bfs(5).options.criterion = 'roc';
%
% options.Xn   = fn;
% options.bcl  = bcl;
% options.bfs  = bfs;
% options.m    = 10;
%
% [bcs,selec] = Bcl_balu(f,d,bcl,bfs,options);
% op.b = bcs; op.v = 10; op.show = 1; op.c = 0.95;
% [p,ci] = Bev_crossval(f(:,selec),d,op);
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [bcs,selec,sp] = Bcl_balu(X,d,bcl,bfs,options)


if nargin == 3
    options = bcl;
    bfs = options.fs;
    bcl = options.cl;
end

if ~isfield(options,'v')
    options.v = 10;
end

if ~isfield(options,'sub')
    options.sub = 1;
end

if ~isfield(options,'outfile')
    options.outfile = 'Bfs_seldata';
end


if ~isfield(options,'c')
    options.c = 0.95;
end

if ~isfield(options,'show')
    options.show = 1;
    % close all
end

if compare(class(bfs),'cell')==0
    bfs = Bfs_build(bfs);
end

if compare(class(bcl),'cell')==0
    bcl = Bcl_build(bcl);
end




f    = X;
fn   = options.Xn;
mmax = options.m;
show = options.show;
v    = options.v;
c    = options.c;

if show
    disp(' ');
    disp(' ');
    disp('Feature and Classifier Selection Tool:')
    disp('======================================')
    % clf
    fsb = Bio_statusbar('Feature selection');
end


f_original = f;
d_original = d;

if options.sub <1
   [f,d] = Bds_stratify(f,d,options.sub);
end

s0      = Bfs_clean(f);
X0      = Bft_norm(f(:,s0),0);
m       = min(mmax,size(X0,2));
n       = length(bcl);   % number of classifiers
nq      = length(bfs);   % number od feature selectors
P       = zeros(n,m,nq);
mf      = 0;             % number of features
qs      = 0;             % number of feature selection algorithm
op.b    = bcl;
op.v    = v;
op.show = 0;
op.c    = c;
sq      = zeros(m,nq);
Sqname  = [];

X = X0;

for q = 1:nq
    if show
        fsb = Bio_statusbar(q/nq,fsb);
    end
    bfs(q).options.m    = m;
    bfs(q).options.show = 0;
    ss = bfs(q).name;
    ft = 0;
    if compare(ss(1:3),'Bfs')
        if compare(ss(1:3),'Bft')==0
            ft = 1;
            bfxname = ss;
        else
            bfxname = ['Bfs_' ss];
        end
    else
        bfxname = ss;

    end
    opf = bfs(q).options;
    if isfield(opf,'b')
        ss = [num2str(q) ':' ss '-' opf.b.name];
    else
        ss = [num2str(q) ':' ss ];
    end
    ss = [ss '                 '];
    ss = ss(1:16);
    if show
        fprintf('\nSelecting %d from %d features (in %d samples) using algorithm %s...\n',m,size(X0,2),size(X0,1),ss)
    end
    if ft==0
        s1 = feval(bfxname,X0,d,bfs(q).options);
    else
        Xs1 = feval(bfxname,X0,d,bfs(q).options);
        Mx = size(X,2);
        Mf = size(f,2);
        Mx1 = size(Xs1,2);
        X = [X Xs1];
        f = [f Xs1];
        s1 = (Mx+1:Mx+Mx1)';
        sf = (Mf+1:Mf+Mx1)';
        s0 = [s0;sf];
        for i=1:Mx1
            si = [bfxname '-' num2str(i) '                             '];
            fn = [fn;si(1:24)];
        end
    end
    fprintf('> %d features selected.\n',length(s1))
    s1 = [s1;zeros(1000,1)];
    Sqname = [Sqname; ss];
    sq (:,q) = s1(1:m);
    save(options.outfile,'Sqname','sq','bfs','options');
end
delete(fsb);

xmf     = [];
ymf     = [];
c1f     = [];
c2f     = [];
STRclass=[];
pkmax   = 0;
if show
    fsb = Bio_statusbar('Classifier selection',fsb);
end
for k=1:m       % features
    if pkmax<1 % when p=100% no more testing!!!

        if show
            fsb = Bio_statusbar(k/m,fsb);
        end

        pqmax = 0;
        if show
            fprintf('\nComputing performances for the %d best feature(s) using...\n',k)
        end
        pnew = 0;
        for q=1:nq  % feature selection algorithms
            s = sq(1:k,q);
            s = s(s>0);
            if length(s)==k
                %            Xk = X(:,s);
                Xk = f_original(:,s0(s));
                if show
                    fprintf('> Feature Selection %s... ',Sqname(q,:))
                end
                [p,ci] = Bev_crossval(Xk,d_original,op);
                P(1:n,k,q) = p;
                for i=1:n
                    cset = bcl(i);
                    % cset.features = s0(s)';
                    cset.ci = ci(i,:);
                    cset.per = p(i);
                    sp(i,k,q) = cset;
                end

                [i,j] = max(p);
                if (i>pqmax)
                    pqmax = i;
                    js = j;
                    qs = q;
                    ssq = s;
                end
                if (pqmax>pkmax)
                    pnew = 1;
                    pkmax = pqmax;
                    c1 = ci(j,1);
                    c2 = ci(j,2);
                    strnclass = [bcl(js).name '            '];
                    strnclass = [strnclass(1:12) ' ' Sqname(qs,:)];
                    mf = k;
                    qss = Sqname(qs,:);
                    ip = js;
                    kp = k;
                    qp = qs;
                    bcs = bcl(js);
                    selec = s0(ssq);
                    st = '*** best ***';
                else
                    st = ' ';
                end
                if show
                    strj = [bcl(j).name '              '];
                    fprintf(' Best classifier %s: Performance = %5.2f%% %s\n',strj(1:8),i*100,st);
                end
            end
        end

        if show
            xmf = [xmf; mf];
            ymf = [ymf;100*pkmax];
            c1f = [c1f;100*c1];
            c2f = [c2f;100*c2];

            plot(xmf,ymf,'b','LineWidth',2);
            plot(xmf,c1f,'r:');
            plot(xmf,c2f,'r:');
            xlabel('Number of features')
            ylabel('Performance [%]');
            axis([0 m 50 100])
            hold on
            pause(0)
            if pnew
                STRclass = [STRclass;sprintf('%3d) %s > %6.2f%%',xmf(end), strnclass,100*max(P(:)))];
            end
        end

    end
end


if show
    delete(fsb);
    legend('Confidence Interval','Performance','Location','SouthEast')
    disp(' ');
    disp(' ');
    disp('Summary:');
    disp('========');
    disp(' ');
    disp('Performance vs. Number of features')
    for i=1:size(STRclass,1)
        fprintf('%s\n',STRclass(i,:))
    end

    disp(' ');
    if ~isempty(fn)
        disp('Selected features for the last classifier:')
        for i=1:length(selec)
            fprintf('%3d) %s\n',i,fn(selec(i),:))
        end
    end
    disp(' ');
    fprintf('                Performance: %5.2f%% in (%5.2f, %5.2f%%) with CI=%5.2f%%\n',P(ip,kp,qp)*100,sp(ip,kp,qp).ci(1)*100,sp(ip,kp,qp).ci(2)*100,c*100);
    fprintf('Folders of Cross Validation: %d\n',v);
    fprintf('                 Classifier: %s\n',bcs.name);
    fprintf('   Feature Selection Method: %s\n',qss);
    fprintf('      %2d Selected Features: %s\n',length(selec),num2str(selec'));
end
bcs.fs_method = qss;