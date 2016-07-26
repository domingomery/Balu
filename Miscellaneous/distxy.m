% function D = distxy(x,y,method)
%
% Toolbox: Balu
%
% It computes distances of each row of x with each row of y
% x is a Nx x n matrix, y is a Ny x n matrix
% D is a Nx x Ny matrix, where D(i,j) = norm(x(i,:)-y(j,:));
% method = 1: it computes row per row
%          2: it computes each row of x with all rows of y
%          3: it computes all rows of x with all rows of y
% method 3 is faster than method 2, and this faster than method 1
% method 3 requires more memory than method 2, and this more than method 1
% method 3 is the default.
%
% Example:
%    x = rand(100,7);
%    y = rand(100,7);
%    D = distxy(x,y); % default is method 3
%
% D.Mery, PUC-DCC, 2012
% http://dmery.ing.puc.cl

function D = distxy(x,y,method)

if ~exist('method','var')
    method = 3;
end


[Nx,nx] = size(x);
[Ny,ny] = size(y);
if nx~=ny
    error('number of elements of columns of x and y must be same...')
end
n = nx;
D = zeros(Nx,Ny);
switch method
    case 1
%        ff = Bio_statusbar('distances');
        for i=1:Nx
%            Bio_statusbar(i/Nx,ff);
            for j=1:Ny
                D(i,j) = norm(x(i,:)-y(j,:));
            end
        end
%        delete(ff)
    case 2
%        ff = Bio_statusbar('distances');
        for i=1:Nx
%            Bio_statusbar(i/Nx,ff);
            di     = ones(Ny,1)*x(i,:)-y;
            D(i,:) = sqrt(sum(di.*di,2)');
        end
%        delete(ff)
    case 3
        xs    = zeros(1,n,Nx);
        xs(:) = x';
        Sx    = repmat(xs,[Ny 1 1]);
        Sx    = Sx-repmat(y,[1 1 Nx]);
        Sx    = Sx.*Sx;
        R     = sqrt(sum(Sx,2));
        Dr    = zeros(Ny,Nx);
        Dr(:) = R;
        D     = Dr';
    case 4
        kd    = vl_kdtreebuild(x');
        [i,d] = vl_kdtreequery(kd,x',y','NumNeighbors',Nx);
        [j,k] = sort(i);
        D = sqrt(d(k));
    case 5 % it is like 3 but with out sqrt
        xs    = zeros(1,n,Nx);
        xs(:) = x';
        Sx    = repmat(xs,[Ny 1 1]);
        Sx    = Sx-repmat(y,[1 1 Nx]);
        Sx    = Sx.*Sx;
        R     = sum(Sx,2);
        Dr    = zeros(Ny,Nx);
        Dr(:) = R;
        D     = Dr';
end

