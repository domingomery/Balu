function [Y,options] = Bft_bow(X,options)


method = lower(options.method(1:2));

switch method
    case {'bu','co'} % build or construct
        fprintf('>>> Building vocabulary with %d words using K-means...\n',options.K);
    case {'vq','freq'} % vector quantification frequencies
        j            = vl_kdtreequery(options.kd,options.Xcen',X','NumNeighbors',1)';
        %[jj, j]      = min(vl_alldist(options.Xcen', options.X'), [], 1) ;

        options.T    = options.Xcen(j,:);
        Y            = (hist(j,1:options.K))';
        h            = Y/sum(Y);
        % options.H    = h;
        options.H    = vl_homkermap(h, 1, 'kchi2', 'gamma', .5) ;
        %Y = options.H;
        

end
