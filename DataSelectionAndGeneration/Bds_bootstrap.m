% [Xb,db,Xnb,dnb] = Bds_bootstrap(X,d,N)
%
% Toolbox: Balu
%
%    Bootstrap sample (with replacement)
%
%    input: (X,d) means features and ideal classification
%    Bbootsample takes a bootstrap sample of (X,d) to build (Xb,db). 
%    N is the number of samples of Xb, default is the number of samples
%    of X.
%    Xnb,dnb are the samples not considered in Xb,db.
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [Xb,db,Xnb,dnb] = Bds_bootstrap(X,d,N)

Nx = size(X,1);
if not(exist('N','var'))
    N = Nx;
end
ib  = mod(ceil(N*rand(N,1)),Nx)+1;
Xb = X(ib,:);
db = d(ib,:);


Nx = size(X,1);
nb = ones(Nx,1);
nb(ib) = 0;
inb = find(nb==1);
Xnb = X(inb,:);
dnb = d(inb,:);

