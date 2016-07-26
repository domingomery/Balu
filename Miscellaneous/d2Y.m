function Y = d2Y(d)
N = length(d);
dmin = min(d);
dmax = max(d);
Y = zeros(N,dmax-dmin+1);
for i=dmin:dmax
    Y(:,i) = d==i;
end
