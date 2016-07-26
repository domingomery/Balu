function p = Bev_confidence(c,e,N)
clf
hist(e)
f = mean(e);
t = (1-c)/2; 
z = norminv(1-t); 
p1 = (f+z^2/2/N-z*sqrt(f/N-f^2/N+z^2/4/N^2))/(1+z^2/N);  
p2 = (f+z^2/2/N+z*sqrt(f/N-f^2/N+z^2/4/N^2))/(1+z^2/N);


z = tinv(1-t,9); 
p1t = (f+z^2/2/N-z*sqrt(f/N-f^2/N+z^2/4/N^2))/(1+z^2/N);  
p2t = (f+z^2/2/N+z*sqrt(f/N-f^2/N+z^2/4/N^2))/(1+z^2/N);



M = length(e);
[i,j] = sort(e);
p1s = e(j(fix(M*t+1.5)));
p2s = e(j(fix(M*(1-t)+0.5)));

p = [p1  f p2
     p1t f p2t
     p1s f p2s];
