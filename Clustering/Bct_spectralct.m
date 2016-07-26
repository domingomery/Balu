% [idx eigvec eigval] = Bct_spectralct(W, k)
%
% Toolbox: Balu
%   Spectral clustering
%
%   idx         -- cluster indexes (same as Matlab k-means).
%   eigvec      -- eigenvectors.
%   eigval      -- eigenvalues.
%
%   W           -- weighted adjacency matrix.
%   k           -- number of clusters.
%
%  Example:
%    load spectraldata
%    d            = ones(size(X,1),1);
%    figure
%    Bio_plotfeatures(X,d)
%    G            = Bct_knngraph2d(X, 50);
%    beta         = 1/0.5^2;
%    [ni nj]      = find(G == true);
%    W            = zeros(size(G));
%    W(G == true) = exp(-beta*(sum((X(ni,:) - X(nj,:)).^2,2)));
%    ds           = Bct_spectralct(W, 3);
%    figure
%    Bio_plotfeatures(X,ds)

%
%
%    See also Bcl_qda.
%
% (c) GRIMA-DCCUC, 2011: Cristobal Moenne
% http://grima.ing.puc.cl

function [idx eigvec eigval] = Bct_spectralct(W, k)
    [eigval eigvec] = eigsmatrix(W, k);

    space = eigvec(:,1:k);
    for ni=1:size(space, 1)
        space(ni,:) = space(ni,:)./norm(space(ni,:));
    end
    idx = kmeans(space, k, 'replicates', 5, 'emptyaction', 'singleton');
end

function [eigval eigvec] = eigsmatrix(W, k)
    D = diag(sum(W));

    [eigvec eigval] = eigs(D\W, k);
    eigval = real(diag(eigval));

    [dummy index] = sort(eigval, 'descend');
    eigval = eigval(index);
    eigvec = eigvec(:,index);
end

