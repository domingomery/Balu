function Ms = Btr_reco3d(kp,P,T)


fra = kp.fra;
img = kp.img;

n = size(P,1)/3;

[N,M] = size(T);

Ms = zeros(4,N);
mm = zeros(3,n);
Pm = zeros(3*n,4);
for i=1:N
    k = 0;
    for j=1:M
        if T(i,j)>0
            k = k+1;
            mm(:,k) = [fra(T(i,j),[2 1]) 1]';
            Pm(indices(k,3),:) = P(indices(img(T(i,j)),3),:);
        end
    end
    ms = mm(:,1:k);
    Ps = Pm(1:3*k,:);
    Mi = Bmv_reco3dn(ms,Ps); %3D affine reconstruction
    Ms(:,i) = Mi;
end
