function d = sqdif(X,Y)
D = X-Y;
D2 = sqrt(sum(D.*D,2));
d = mean(D2);