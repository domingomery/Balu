function bfs = Bfs_build(varargin)

v = varargin;
n = nargin;
if compare(class(v),'cell')==0
    v = v{1};
    n = length(v);
end

for k=1:n
    s = lower(char(v(k)));
    ns = length(s);
    ok = 0;
    if ns>7
        if compare(s(1:7),'sfs-knn')==0
            bfs(k).name = 'sfs';
            bfs(k).options.b.name    = 'knn';
            bfs(k).options.b.options.k = str2num(s(8:end));
            ok = 1;
        end
        if compare(s(1:7),'sfs-svm')==0
            bfs(k).name = 'sfs';
            bfs(k).options.b.name    = 'svmplus';
            bfs(k).options.b.options.kernel = str2num(s(8:end));
            ok = 1;
        end
        if ns>9
            if compare(s(1:9),'sfs-nnglm')==0
                bfs(k).name = 'sfs';
                bfs(k).options.b.name    = 'nnglm';
                bfs(k).options.b.options.method = str2num(s(10:end));
                bfs(k).options.b.options.iter = 12;
                ok = 1;
            end
        end
    end

    if not(ok)
        switch s
            case 'fosmod';
                bfs(k).name = 'fosmod';
                bfs(k).options = [];

            case 'lsef';
                bfs(k).name = 'lsef';
                bfs(k).options = [];

            case 'sfs-fisher';
                bfs(k).name = 'sfs';
                bfs(k).options.b.name    = 'fisher';

            case 'sfs-sp100';
                bfs(k).name = 'sfs';
                bfs(k).options.b.name    = 'sp100';

            case 'sfs-lda'
                bfs(k).name = 'sfs';
                bfs(k).options.b.name    = 'lda';
                bfs(k).options.b.options.p = [];
            case'sfs-qda'
                bfs(k).name = 'sfs';
                bfs(k).options.b.name    = 'qda';
                bfs(k).options.b.options.p = [];
            case 'rank-roc';
                bfs(k).name = 'rank';
                bfs(k).options.criterion = 'roc';

            case 'rank-ttest';
                bfs(k).name = 'rank';
                bfs(k).options.criterion = 'ttest';

            case 'rank-entropy';
                bfs(k).name = 'rank';
                bfs(k).options.criterion = 'entropy';

            case 'plsr';
                bfs(k).name = 'Bft_plsr';

            case 'pca';
                bfs(k).name = 'Bft_pca';

            case 'lsef-plsr';
                bfs(k).name = 'Bft_lseft';
                bfs(k).options.pca = 0;

            case 'lsef-pca';
                bfs(k).name = 'Bft_lseft';
                bfs(k).options.pca = 1;

            case 'mrmr';
                bfs(k).name = 'Bfs_mRMR';

            case 'all';
                bfs(k).name = 'Bfs_all';

            otherwise
                error('Bfs_build does not recognize %s as feature selection method.',s)

        end
    end
end