function R = Bmv_line2img(ell,NM)
N = NM(1);
M = NM(2);
R = zeros(N,M);
if ell(1)~=0
    for j=1:M
        i = round(-(ell(2)*j+ell(3))/ell(1));
        if and(i>0,i<=N)
            R(i,j) = 1;
        end
    end
end
if ell(2)~=0
    for i=1:N
        j = round(-(ell(1)*i+ell(3))/ell(2));
        if and(j>0,j<=M)
            R(i,j) = 1;
        end
    end
end