% selec = Bfs_random(X,d,options)
%
% Toolbox: Balu
%    Select best features from random subsets of X according to ideal
%    classification d. 
%    options.m is the number features will be selected.
%    options.M is the number of random subsets to be tested.
%    options.b.method = 'fisher' uses Fisher objetctive function.
%    options.b.method = 'sp100' uses as criteria Sp @Sn=100%.
%    options.b can be any Balu classifier structure (see example).
%    options.show = 1 display results.
%    selec is the indices of the selected features.
%
% Example 1: Feature selection using Fisher dicriminant
%    load datareal
%    op.m = 10;                     % 10 features will be selected
%    op.M = 1000;                   % number of random sets to be tested
%    op.show = 1;                   % display results
%    op.b.name = 'fisher';          % Feature score
%    s = Bfs_random(f,d,op);        % index of selected features
%    X = f(:,s);                    % selected features
%    Xn = fn(s,:)                   % list of feature names
%
% Example 2: Feature selection using a KNN classifier
%    load datareal
%    op.m = 10;                     % 10 features will be selected
%    op.M = 1000;                   % number of random sets to be tested
%    op.show = 1;                   % display results
%    op.b.name = 'knn';             % Feature score is the performance
%    op.b.options.k = 5;            % 5 neighbors
%    s = Bfs_random(f,d,op);        % index of selected features
%    X = f(:,s);                    % selected features
%    Xn = fn(s,:)                   % list of feature names
%
% (c) D.Mery, PUC-DCC, Jul. 2011
% http://dmery.ing.puc.cl

function selec = Bfs_random(X,d,options)


m = options.m;
M = options.M;
show = options.show;

s0 = Bfs_clean(X);

f = X(:,s0);
nf = size(f,2);
Jmax = 0;

k = 0;
while (k<=M)
    N = rand(nf,1);
    [i,j] = sort(N);
    s = j(1:m);
    fs = f(:,s);
    Js = Bfa_score(fs,d,options);
    if (Js>Jmax)
        selec = s;
        Jmax = Js;
        if show
            fprintf('Jmax = %8.4f\n',Jmax);
        end
    end
    k = k+1;
end
