function Y = Bim_maxmin(X)
X = double(X);
Xmax = max2(X);
Xmin = min2(X);
Y = (X-Xmin)/(Xmax-Xmin);