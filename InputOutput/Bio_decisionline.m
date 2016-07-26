% Bio_decisionline(X,d,op)
%
% Toolbox: Balu
%
%    Diaplay a 2D feature space and decision line.
%
%    X: Sample data
%    d: classification of samples
%    op: output of a trained classifier.
%
%    Example:
%       load datagauss                             % simulated data (2 classes, 2 features)
%       Xn = ['\beta_1';'\beta_2'];
%       b(1).name = 'knn';   b(1).options.k = 5;                              % KNN with 5 neighbors
%       b(2).name = 'lda';   b(2).options.p = [];                             % LDA
%       b(3).name = 'svm';   b(3).options.kernel = 4;                         % rbf-SVM
%       op = b;
%       op = Bcl_structure(X,d,op);
%       close all
%       Bio_decisionline(X,d,Xn,op);
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function Bio_decisionline(X,d,Xn,op)
if size(X,2)~=2
    error('Bio_decisionline works for two features only.')
end
clf
Bio_plotfeatures(X,d,Xn);
n = length(op);
hold on
ax = axis;
s=0.1;
x = (ax(1)):s:(ax(2));
nx = length(x);
y = (ax(3)):s:(ax(4));
ny = length(y);
rx = ones(ny,1)*x;
ry = y'*ones(1,nx);
XXt = [rx(:) ry(:)];
op = Bcl_structure(X,d,op);
dt = Bcl_structure(XXt,op);

dmin = min(d);
dmax = max(d);
scol = 'ycwkbr';
tcol = 'brgywk';
close all
for k=1:n
    figure
    Bio_plotfeatures(X,d,Xn);
    title(op(k).options.string)
    for i=dmin:dmax
        ii = find(dt(:,k)==i);
        plot(XXt(ii,1),XXt(ii,2),[scol(i-dmin+1) 'o']);
        ii = find(d==i);
        plot(X(ii,1),X(ii,2),[tcol(i-dmin+1) 'x']);
    end
end

