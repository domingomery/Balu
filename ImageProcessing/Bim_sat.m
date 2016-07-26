function y = Bim_sat(x,xmin,xmax)
y = x;
imin = find(x<xmin);
imax = find(x>xmax);
if not(isempty(imin))
    y(imin) = xmin;
end
if not(isempty(imax))
    y(imax) = xmax;
end
