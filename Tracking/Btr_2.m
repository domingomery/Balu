% B = Btr_2(kp,F,options)
%
% Toolbox: Balu
%
%    Matching points between all two views p and q of a sequence, for
%    p=1:n-1, and for q=p+1:p+m (n is the number of views in the sequence
%    and m is defined by options.mviews)
%
%    kp keypoints structure according function Bsq_des (see help)
%
%    F are the fundamental matrices according Bsq_fundamental (see help).
%
%    B is a Nx2 matrix with N matchings. a = B(i,1) and b = B(i,2) mean that
%    a matching between keypoints a and b was found. The position of these
%    keypoints are ma = kp.fra(a,[2 1]) and mb = kp.fra(b,[2 1]) in views
%    p = kp.img(a) and q = kp.img(b). The descriptions of these keypoints 
%    (feature vectors, e.g., SIFT) are xa = kp.des(a) and xb = kp.des(b).
%
%    A matching is found if the following constraints are fullfil:
%       1) epipolar distance (ma,mb) < options.depimax
%       2) Xi^2 distance (xa,xb)     < options.dxi2max
%       3) norm(ma-mb)/(q-p)         < options.dnmax
%
%
%  Example:
%      See example in Btr_demo.
%
%  See also Bsq_des, Bsq_fundamental, Btr_3, Btr_demo.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function B = Btr_2(kp,F,options)


dxi2max = options.dxi2max;
dnmax   = options.dnmax;
depimax = options.depimax;
m       = options.mviews;

img     = kp.img;
des     = kp.des;
fra     = kp.fra;

p1      = min(kp.img);
p2      = max(kp.img);

B = zeros(20000,2);
ib = 0;
N = max(img);
e = zeros(N,1);
f = zeros(N,1);

for i=1:N
    ii = find(img==i);
    e(i) = min(ii);
    f(i) = max(ii);
end

Fpq = zeros(3,3);
for p=p1:p2
    for a = e(p):f(p)
        xa = double(des(a,:));
        ma = fra(a,[2 1]);
        for q = p+1:min([p+m N])
            Fpq(:) = F(p,q,:);
            for b = e(q):f(q)
                mb   = fra(b,[2 1]);
                d    = norm(ma-mb)/(q-p);
                if (d<dnmax)
                    xb   = double(des(b,:));
                    dXi  = Bfa_dXi2(xa,xb);
                    if (dXi<dxi2max)
                        dEpi   = Bmv_epidist([ma 1]',[mb 1]',Fpq);
                        if (dEpi<depimax)
                            ib      = ib+1;
                            B(ib,:) = [a b];
                        end
                    end
                end
            end
        end
        
    end
end
B = B(1:ib,:);

if options.show
    fprintf('Btr_2      : %4d matchings in 2 views.\n',size(B,1))
end
    
