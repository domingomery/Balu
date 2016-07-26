function I = Bsq_patch(f,images,m,d,r)
n = length(images);
if size(m,1)==3
    m = h2i(m);
end
I = zeros(r,n*r);
for i=1:n
    Ii = Bio_loadimg(f,images(i));
    [N,M] = size(Ii);
    Ij = zeros(N+2*d,M+2*d);
    Ij(d+1:d+N,d+1:d+M) = Ii;
    mi = round(m(1,i)-d:m(1,i)+d)+d;
    mj = round(m(2,i)-d:m(2,i)+d)+d;
    Pj = imresize(Ij(mi,mj),[r r]);
    I(:,indices(i,r)) = Pj;
end
