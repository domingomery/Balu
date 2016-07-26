% selec = Bfs_balu(X,d,options)
%
% Toolbox: Balu
%    Feature selection of "best" options.m features of X to ideal
%    classification d. This function uses only a portion of options.s
%    samples to select the features. 
%    Bfs_balu (1) normalizes (using Bft_norm), (2) cleans (using Bfs_clean)
%    and (3) selects features (using Bfs_sfs).
%    options.b.method = 'fisher' uses Fisher objetctive function.
%    options.b.method = 'sp100' uses as criteria Sp @Sn=100%.
%    options.b can be any Balu classifier structure (see example).
%    options.show = 1 display results.
%    selec is the indices of the selected features.
%
% Example 1: Balu feature selection using Fisher dicriminant
%    load datareal
%    op.m = 10;                     % 10 features will be selected
%    op.s = 0.75;                   % only 75% of sample will be used
%    op.show = 1;                   % display results
%    op.b.name = 'fisher';          % SFS with Fisher
%    s = Bfs_balu(f,d,op);          % index of selected features
%    X = f(:,s);                    % selected features
%    Xn = fn(s,:)                   % list of feature names
%
% Example 2: Balu feature selection using a KNN classifier
%    load datareal
%    op.m = 10;                     % 10 features will be selected
%    op.s = 0.75;                   % only 75% of sample will be used
%    op.show = 1;                   % display results
%    op.b.name = 'knn';             % SFS with KNN
%    op.b.options.k = 5;            % 5 neighbors
%    s = Bfs_balu(f,d,op);          % index of selected features
%    X = f(:,s);                    % selected features
%    Xn = fn(s,:)                   % list of feature names
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function selec = Bfs_balu(X,d,options)

s = options.s;
show = options.show;
m = options.m;


if (s<0)||(s>1)
    error('Bfs_balu: options.s must be between 0 and 1...')
end

fT = X;
dT = d;

if show
    close all
end

if s<1
    [f,d] = Bds_stratify(fT,dT,s);
else
    f = X;
end


% cleaning
if show
    disp('eliminating no relevant and correlated features...');
end
selec0 = Bfs_clean(f);
f = f(:,selec0);

% normalize
if show
    disp('normalizing features...');
end
ff = Bft_norm(f,1);

% selection
if show
    disp('selecting features...');
end

M = size(ff,2);
if M<m
    error('Bfs_balu: number of features to be selected (%d) is larger than number of existing features (%d)',M,m);
end
warning off %#ok<WNOFF>
N = nchoosek(M,m);
warning on %#ok<WNON>

if N>10000
    if show
        disp('selecting features using SFS...');
    end
    selec = Bfs_sfs(ff,d,options);
else
    if show
        disp('selecting features using exhaustive search...');
    end
    selec = Bfs_exsearch(ff,d,options);
end

selec = selec0(selec);
