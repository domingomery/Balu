% [P,H1] = Btr_sfm(kp,Ho,options)
%
% Toolbox: Balu
%
%    Structure from Motion.
%
%    kp keypoints structure according function Bsq_des (see help)
%
%    Ho is a matching multi-views matrix with Nxn indices for N matchings
%    in n views.
%
%    options.sfs_method = 1 for affine projection and 2 for projective
%    reconstruction.
%
%    options.sfm_samples is the number of samples to be considered by
%    adjustment (minimum is 4). If sfm_samples is 0 then all samples will
%    taken into account and no RANSAC will be performed.
%
%    options.sfm_iter is the number of RANSAC iterations, for sfs_samples=0
%    use sfm_iter = 1.
%
%    The inliers must be cover at least options.sfm_cover in each dimension
%    of the image (e.g. 0.5).
%
%    P includes the projection matrices of n views as follows:
%    Projection Pk = P(k*3-2:k*3,:), for k=1,...,n
%
%    H1 is the matching multi-view matrix with the inliers matchings.
%
%  Example:
%      See example in Btr_demo.
%
%  See also Bsq_des, Btr_join, Btr_siftn, Btr_demo.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [P,H1] = Btr_sfm(kp,Ho,options)


show       = options.show;

nq = options.sfm_samples;
cv = options.sfm_cover;

x   = kp.fra(Ho,1);
y   = kp.fra(Ho,2);
dxx = (max(x)-min(x))*cv;
dyy = (max(y)-min(y))*cv;

ok  = 0;
k   = 0;

while not(ok)
    k = k + 1;
    [P,H1] = Btr_ba(kp,Ho,options); % function included below in this file
    x  = kp.fra(H1,1);
    y  = kp.fra(H1,2);
    dx = max(x)-min(x);
    dy = max(y)-min(y);
    if (nq == 0) || (cv == 0)
        ok = 1;
    else
        % the inliers must cover at least 50% of each dimension of the
        % image
        ok = (size(H1,1)>=4) & (dx>dxx) & (dy>dyy);
    end
    if show
        fprintf('Btr_sfm    : %4d trajectories found, dx = %f (%f), dy=%f (%f).\n',size(H1,1),dx,dxx,dy,dyy);
    end
end

end

function [Ps,H1] = Btr_ba(kp,Ho,options)


method = options.sfm_method;
Nq     = options.sfm_samples;
Nk     = options.sfm_iter;

img = kp.img;
fra = kp.fra;


m = size(Ho,2);
n = size(Ho,1);

g = min(kp.img):max(kp.img);

R = sum(ones(n,1)*g-img(Ho),2);

Ho = Ho(R==0,:);
n = size(Ho,1);

if Nq==0
    Nq = n;
end

if n>=4
    xt = zeros(2,n,m);
    for i=1:m
        for j=1:n
            xt(:,j,i) = fra(Ho(j,i),[2 1])';
        end
    end

    % RANSAC
    % i = 1...m views
    % j = 1...n 3D points
    % q random samples

    q       = Nq;  % random samples
    emin    = 5;   % minimum distance for inlieres
    outlmin = Inf;
    k       = 1;
    while (k<=Nk) && (outlmin>0)
        j = vl_colsubset(1:n,q);
        x = xt(:,j,:);
        if method==1
            [Xs,P] = Bmv_bundleafin(x);
            inl  = zeros(n,1);
            for j=1:n
                mm = ones(3,m);
                for i=1:m
                    mm(1:2,i) = xt(:,j,i);
                end
                [Mj,e] = Bmv_reco3dna(mm,P); %3D affine reconstruction
                inl(j) = max(e)<emin;
            end
        else
            [n1,n2,n3] = size(x);
            xx = ones(n1+1,n2,n3);
            xx(1:2,:,:) = x;
            [Xs,P] = Bmv_bundleproj(xx);
            inl  = zeros(n,1);
            for j=1:n
                mm = ones(3,m);
                for i=1:m
                    mm(1:2,i) = xt(:,j,i);
                end
                [Mj,e] = Bmv_reco3dn(mm,P); %3D projective reconstruction
                inl(j) = max(e)<emin;
            end
        end
        outl = n-sum(inl);
        if outl<outlmin
            outlmin = outl;
            inlmin = inl;
            Ps  = P;
            % fprintf('%d/%d outliers = %d/%d\n',k,Nk,outl,n)
        end
        k = k+1;
    end
    H1 = Ho(inlmin==1,:);
else
    H1 = -1;
    Ps = -1;
    fprintf('Error: bundle adjustment requieres 4 or more points, there is only %d.\n',n);
end
end