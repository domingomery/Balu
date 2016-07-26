function d = Bds_labels(v)
d = zeros(sum(v),1,'uint8');
k = 0;
for i=1:length(v)
    d(k+1:k+v(i))=i;
    k = k+v(i);
end