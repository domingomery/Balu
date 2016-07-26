% selec = Bfs_exsearch(X,d,options)
%
% Toolbox: Balu
%    Feature selection using exhaustive search for fatures X according to 
%    ideal classification d. optins.m features will be selected.
%    options.b.method = 'fisher' uses Fisher objetctive function.
%    options.b.method = 'sp100' uses as criteria Sp @Sn=100%.
%    options.b can be any Balu classifier structure (see example).
%    options.show = 1 display results.
%    selec is the indices of the selected features.
%
% Example 1: Exhaustive search from using Fisher dicriminant
%    load datareal
%    s1 = [279 235 268 230 175 165 207 160 269 157]; %indices using Example
%                                                    %of Bfs_sfs.
%    X1 = f(:,s1);                    % preselected features
%    Xn1 = fn(s1,:);
%    op.m = 3;                        % 3 features will be selected
%    op.show = 1;                     % display results
%    op.b.name = 'fisher';            % SFS with Fisher
%    s = Bfs_exsearch(X1,d,op);       % index of selected features
%    X2 = X1(:,s);                    % selected features
%    Xn2 = Xn1(s,:)                   % list of feature names
%
% Example 2: Exhaustive serach using a KNN classifier
%    load datareal
%    s1 = [279 235 268 230 175 165 207 160 269 157]; %indices using Example
%                                                    %of Bfs_sfs.
%    X1 = f(:,s1);                    % preselected features
%    Xn1 = fn(s1,:);
%    op.m = 4;                        % 4 features will be selected
%    op.show = 1;                     % display results
%    op.b.name = 'knn';               % Feature selection using KNN
%    op.b.options.k = 5;              % 5 neighbors
%    s = Bfs_exsearch(X1,d,op);       % index of selected features
%    X2 = X1(:,s);                    % selected features
%    Xn2 = Xn1(s,:)                   % list of feature names
%
% D.Mery, PUC-DCC, Jul. 2011
% http://dmery.ing.puc.cl

function selec = Bfs_exsearch(X,d,options)

m = options.m;
show = options.show;
f = X;

M = size(f,2);

N = nchoosek(M,m);

if (N>10000)
    ok = input(sprintf('Exhaustive Search needs %d evaluations... continue [yes=1/no=0]?',N));
    if not(ok)
        error('Exhaustive search for feature selection interrupted.')
    end
end

T = nchoosek(1:M,m);

Jmax = 0;
for i=1:N
    fs = f(:,T(i,:));
    Js = Bfa_score(fs,d,options);
    if (Js>Jmax)
        Jmax = Js;
        ks = i;
        if show
           fprintf('step=%2d/%d J=%f\n',i,N,Jmax)
        end
    end
end
selec = T(ks,:)';
