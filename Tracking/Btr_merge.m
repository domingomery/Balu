% E = Btr_merge(kp,D)
%
% Toolbox: Balu
%
%    Merge tracks with common matching points.
%
%    kp keypoints structure according function Bsq_des (see help)
%
%    D is a Nxm matrix with N matchings in m views. The output E is a
%    matrix with merged trajectories with common keypoints.
%
%  Example:
%      See example in Btr_demo.
%
%  See also Bsq_des, Bsq_fundamental, Btr_2, Btr_3, Btr_join, Btr_demo.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function E = Btr_merge(kp,D)

m = max(kp.img);

T = zeros(1000,m);
d = D(1,:);
d = d(d>0);
nd = length(d);
T(1,1:nd) = d;
n=1;
N4 = size(D,1);
for i=2:N4
    tn = 0;
    d = D(i,:);
    d = d(d>0);
    nd = length(d);
    for j=1:n
        sd = 0;
        for k=1:nd
            sd = sum(d(k)==T(j,:))+sd;
        end
        if sd>=(nd-1)
            s = sort(unique([T(j,:) D(i,:)]));
            s = s(s>0);
            [si,sj] = unique(kp.img(s));
            s = s(sj);
            ns = length(s);
            T(j,1:ns) = s;
            tn = 1;
        end
    end
    if tn==0;
        n=n+1;
        s = D(i,:);
        s = s(s>0);
        T(n,1:length(s)) = s;
    end
end
T = T(1:n,:);
s = sum(T,1);
T = T(:,s>0);



e = ones(n,1);
for i=1:n
    H = ones(n,1)*T(i,:)-T;
    H2 = H.*H;
    D2 = sum(H2,2);
    ii = find(abs(D2)<1e-5);
    if length(ii)>1
        ii(1) = [];
        e(ii) = 0;
    end
end
E = T(e==1,:);
