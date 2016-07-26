% C = Btr_3(kp,B,T,options)
%
% Toolbox: Balu
%
%    Matching points between all three views p, q and r of a sequence, for
%    p=1:n-1, and for q=p+1:p+m (n is the number of views in the sequence
%    and m is defined by options.mviews)
%
%    kp keypoints structure according function Bsq_des (see help)
%
%    B is a matching two views matrix according function Btr_2 (see help)
%
%    T are the trifocal tensors according Bsq_trifocal (see help).
%
%    C is a Nx3 matrix with N matchings. a=C(i,1),  b=B(i,2) and c=B(i,3)
%    mean that a matching between keypoints a, b and cwas found. The position
%    of these keypoints are ma=kp.fra(a,[2 1]), mb=kp.fra(b,[2 1]) and
%    mc=kp.fra(c,[2 1]) in views p=kp.img(a), q=kp.img(b) and r=kp.img(c).
%
%    A matching is found if the reproyection of mc estimated using ma and
%    mb (and the trifocal tensors) is smaller than option.dtrimax.
%
%
%  Example:
%      See example in Btr_demo.
%
%  See also Bsq_des, Bsq_fundamental, Btr_2, Btr_demo.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function C = Btr_3(kp,B,T,options)

dtrimax2 = (options.dtrimax)^2;

options.join_iter = 1;
C = Btr_join(B,[],options);

if ~isempty(C)
    mg = max(kp.img);
    Ts = zeros(mg^3,27);
    mg2 = mg^2;
    mg0 = mg2+mg;
    for p=1:mg-2
        for q=p+1:mg-1
            for r=q+1:mg
                Tpqr = T(p,q,r,:);
                Ts(p*mg2+q*mg+r-mg0,:) = Tpqr;
            end
        end
    end
    
    if dtrimax2 < 1e6
        mg = max(kp.img);
        mg2 = mg^2;
        mg0 = mg2+mg;
        
        nc = size(C,1);
        m0 = zeros(3*nc,3);
        m0(:) = [double(kp.fra(C(:,:),[2 1])) ones(3*nc,1)]';
        
        nt    = [kp.img(C) ones(nc,1)]*[mg2 mg 1 -mg0]';
        
        m1    = zeros(3,nc); m1(:) = m0(:,1);
        m2    = zeros(3,nc); m2(:) = m0(:,2); m21 = (ones(3,1)*m2(1,:))';
        m3    = zeros(3,nc); m3(:) = m0(:,3);
        
        m3s   = [sum((Ts(nt, 1: 3) -m21.*Ts(nt, 7: 9))'.*m1)
            sum((Ts(nt,10:12) -m21.*Ts(nt,16:18))'.*m1)
            sum((Ts(nt,19:21) -m21.*Ts(nt,25:27))'.*m1)];
        
        m3s   = m3s./(ones(3,1)*m3s(3,:));
        
        d     = m3(1:2,:)-m3s(1:2,:);
        d2    = (sum(d.*d));
        
        C = C(d2<dtrimax2,:);
    end
end
if options.show
    fprintf('Btr_3      : %4d matchings in 3 views.\n',size(C,1))
end






