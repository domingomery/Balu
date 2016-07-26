% H2 = Btr_join(H0,H1,options)
% H2 = Btr_join(H0,[],options)
% H2 = Btr_join(H0,p,options)
%
% Toolbox: Balu
%
%    Join of matching points.
%
%    Hk (k=0,1,2) is a matching multi-views matrix with Nkxnk indices for
%    Nk matchings in nk views. H2 has the tracks that have p common
%    elements in H1 and H0 (the last p elements of a row in H0 must be
%    equal to the first p elements of a row in H1). p is defined as n1-1
%    (if H1 is given). If H1 is empty H1=H0.
%
%    options.show display results
%    options.join_iter is the number of repetitions of this process, for
%    values grater than 1, H0 = Btr_join(H0,H1,options) will be repeated,
%    finally H2 is the resulting H0.
%
%  Example:
%      See example in Btr_demo.
%
%  See also Bsq_des, Btr_join, Btr_siftn, Btr_demo.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


function H2 = Btr_join(H0,H1,options)


if isempty(H0)
    H2 = [];
else
    
    
    nj = options.join_iter;
    
    if nj>1
        options.join_iter = 1;
        for i=1:nj
            H0 = Btr_join(H0,H1,options);          % Structure Tracks
        end
        H2 = H0;
    else
        
        
        m = size(H0,2);
        
        if isempty(H1)
            H1 = m-1;
        end
        
        if length(H1(:))==1
            p   = H1;
            X   = H0(:,1:p);
            Xt  = H0(:,m-p+1:m);
            Ho  = H0;
        else
            p   = size(H1,2)-1;
            X   = H1(:,1:p);
            Xt  = H0(:,m-p+1:m);
            Ho  = H1;
        end
        mo = size(Ho,2);
        
        
        nd = 30;
        kdtree = vl_kdtreebuild(X');
        [q,d]  = vl_kdtreequery(kdtree,X',Xt','NumNeighbors',nd);
        
        [ii,jj] = find(d==0);
        kk = jj*nd+ii-nd;
        
        nk = length(ii);
        
        if nk>0
            i = jj;j=q(kk);
            H2 = [H0(i,:) Ho(j,p+1:mo)];
        else
            H2 =[];
        end
        if options.show
            fprintf('Btr_join   : %4d matchings in %d views.\n',nk,size(H2,2))
        end
    end
end