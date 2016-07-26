function G = Bct_knngraph2d(data,k)

n = size(data,1);
G = false(n,n);

for ni=1:n
    d = ones(n,1)*data(ni,:) - data;
    dd = d.*d;
    e = sum(dd,2);
    
    [a b] = sort(e);
    v = b(2:k+1);
    for nj=1:k
        G(ni,v(nj)) = true;
    end
end
G = logical(~diag(ones(n,1)).*G + diag(ones(n,1)));
end