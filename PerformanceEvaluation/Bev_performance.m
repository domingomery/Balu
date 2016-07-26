% p = Bev_erformance(d1,d2,nn)
%
% Toolbox: Balu
%    Performance evaluation between two classifications, e.g., ideal (d1) 
%    and real (d2) classification.
%
%    d1 and d2 are vectors or matrices (vector Nxn1 and Nxn2 respectivelly)
%    at least n1 or n2 muts be one. N is the number of samples.
%    p is the performance, diagonal sum divided by N.
%    The classes should be labeled as 1, 2, ... n, but
%    if nn = [i j] is given it indicates that the lowest
%    class is min(nn) and the highest class is max(nn).
%    Thus, n = max(nn)-min(nn)+1 and the size of T is nxn.
%
%    d1 (or d2) can be a matrix (Nxn1 or Nxn2) containing the classification 
%    of n1 (or n2) classifiers. In this case the output p is a vector 
%    n1x1 (or n2x1) containing the performance of each classifier. 
%    If n1=1 and n2=1, then p is a scalar.
%
%    Bperformance(d1,d2) and Bperformance(d2,d1) obtain the same result.
%
%    Example:    
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.p = [];
%       ds1 = Bcl_lda(X,d,Xt,op);  % LDA classifier
%       ds2 = Bcl_qda(X,d,Xt,op);  % QDA classifier
%       ds = [ds1 ds2];
%       p = Bev_performance(ds,dt) % performance on test data
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function p = Bev_performance(d1,d2,nn)

d1 = single(d1);
d2 = single(d2);

n1 = size(d1,2);
n2 = size(d2,2);

if (n1==1) || (n2==1)

    if n2>n1
        ds = d2;
        d  = d1;
    else
        d  = d2;
        ds = d1;
    end

    n = size(ds,2);
    p = zeros(n,1);
    for i=1:n
        if exist('nn','var');
            [T,pp] = Bev_confusion(d,ds(:,i),nn);
        else
            [T,pp] = Bev_confusion(d,ds(:,i));
        end
        p(i) = mean(pp);
    end

else
    error('Bev_performance: at least d1 or d2 must have only one column.');
end