function d = compare(X,Y)

%[Nx,Mx] = size(X);
%[Ny,My] = size(Y);

%if or(Nx~=Ny,Mx~=My)
%    fprintf('Warning: sizes are not equal (%d,%d) and (%d,%d).\n',Nx,Mx,Ny,My)
%end

x = double(X(:));
y = double(Y(:));

nx = length(x);
ny = length(y);

nd = abs(nx-ny);

if nd==0
    D = x-y;
    d = sum(sqrt(D.*D));
else
    d = -nd;
end