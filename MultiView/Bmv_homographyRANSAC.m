% H = Bmv_homographyRANSAC(m1,m2)
%
% Toolbox: Balu
%
%    Estimation of Homography Matrix using RANSAC.
%
%    m1 and m2 are n corresponding points in two views (m1(:,k) and m2(:,k)
%    are the k-th corresponding points (for k=1..n) stored as homogeneous
%    coordinates.
%
%    Example:
%       m1 = [rand(2,20);ones(1,20)];           % projection points in view 1 
%       H = rand(3,3);                          % original homography matrix
%       H = H/H(3,3)                            % normalization
%       m2 = H*m1;m2 = m2./(ones(3,1)*m2(3,:)); % projection points in view 2
%       m2(:,20) = [10 10 1]';                  % outlier
%       Hs = Bmv_homographyRANSAC(m1,m2);       % estimations of H
%       Hs = Hs/Hs(3,3)                         % normalization
%       m2s = Hs*m1;m2s = m2s./(ones(3,1)*m2s(3,:));
%       e = m2s-m2                              % reprojection error
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function H = Bmv_homographyRANSAC(ma,mb)

n = size(ma,2);

np = 4;

xa = ma(1:2,:);
xb = mb(1:2,:);

t = max([xa(:);xb(:)])/20;

if n>=np

    our = Inf;
    for k=1:15000
        qq = rand(n,1);
        [ii,jj] = sort(qq); %#ok<ASGLU>
        j = jj(1:np);
        m1 = ma(:,j);
        m2 = mb(:,j);
        H = Bmv_homographySVD(m1,m2);

        mbs = H*ma; mbs = mbs./(ones(3,1)*mbs(3,:));
        mas = H\mb; mas = mas./(ones(3,1)*mas(3,:));
        d1 = mbs(1:2,:)'-mb(1:2,:)';
        d2 = mas(1:2,:)'-ma(1:2,:)';
        dm = sqrt(sum(d1.*d1,2))+sqrt(sum(d2.*d2,2));
        ii = find(dm>t);
        o = length(ii);

        if (o<our)
            our=o;
            Hmin = H;
        end

    end
    H = Hmin;
else
    error('Bmv_homographyRANSAC: not enough matching points');
end