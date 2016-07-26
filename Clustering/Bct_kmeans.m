% function [ds,C] = Bct_kmeans(X,k,show)
%
% Toolbox: Balu
%
%    k-means clustering
%    X matrix of samples
%    k number of clusters
%    ds assigned class number
%    C centroids of each cluster
%    show = 1 display intermediate results
%    show = 0 uses kmeans algortithm of vlfeat (if installed else algorithm
%    of matlab).
%
%    Example:
%     [X,d] = Bds_gaussgen([10 1;1 10;15 15],4*ones(3,3),100*ones(3,1));
%     figure(1)
%     Bio_plotfeatures(X,d);
%     ds = Bct_kmeans(X,3);
%     figure(2)
%     Bio_plotfeatures(X,ds);
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [ds,C] = Bct_kmeans(X,k,show)

if not(exist('show','var'))
    show = 0;
end


X = double(X);

if show==0
    if ~exist('vl_kmeans','file')

        [ds,C] = kmeans(X,k);
    else
        % [ds,C] = vl_kmeans(X,k);
        [C,ds] = vl_kmeans(X',k);
        ds = ds'; C = C';
    end
else
    fprintf('Computing K-means clustering for %d clusters and %dx%d data...\n',k,size(X,1),size(X,2));


    [N,P] = size(X);
    d = fix(k*0.999*rand(N,1)+1);
    ok = 0;

    C = zeros(k,P);
    while not(ok)
        Dk = Inf*ones(N,k);
        for i=1:k
            ii = find(d==i);
            if isempty(ii)
                mi = rand(1,P);
            else
                mi = mean(X(ii,:),1);
            end
            D       = X-ones(N,1)*mi;
            D2      = D.*D;
            Dk(:,i) = sum(D2,2);
            C(i,:)  = mi;
        end
        [i,j] = min(Dk,[],2);
        ds = j;

        e = norm(d-ds);
        d = ds;
        if e<1
            ok = 1;
        end
        if show
            if P<=3
                clf
                Bio_plotfeatures(X,ds)
            else
                clf
                Bio_plotfeatures(X(:,1:3),ds)
            end
            pause(1)
        end
    end
end
