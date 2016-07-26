% [Z,f] = Btr_analysis(W,Iw,options)
%
%
% Toolbox: Balu
%
%    Track analysis. Btr_analysis selects those trajectories that
%
%    The 3D reconstructed points are reprojected in those views where the
%    segmentation may have failed to obtain the complete track in all
%    views. The reprojected points should correspond to the centroids of
%    the non-segmented regions. We calculate the size of the projected
%    region as an average of the sizes of the identified regions in the
%    track. In each view, we define a small window centered in the computed
%    centroids. Afterwards, we average the tracked windows and the tracked
%    Harris windows. Since regions must appear as contrasted zones relating
%    to their environment, we verify if the contrast of each averaged
%    window is large enough. The windows can be rectified using a set of
%    homographies which represent possible local image transformation.
%
%    kp keypoints structure according function Bsq_des (see help)
%
%    Y is a Nxm matrix with N matchings in m views.
%
%    P includes the projection matrices of n views as follows:
%    Projection Pk = P(k*3-2:k*3,:), for k=1,...,m
%
%    files is a structure that define the images of the sequence according
%    to function Bloadimg (see help).
%
%
%  Example:
%      See example in Btr_demo.
%
%  See also Bsq_des, Btr_join, Btr_merge, Btr_2, Btr_3, Btr_demo.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [Z,Iz,Wd] = Btr_classify(W,Iw,Y,options)


n = size(W,3);
d = ones(n,1);
m = size(W,1);
p = size(Iw,1)/n;
ny = size(Y,2);

Iz = Iw;
nz = ones(size(Iz,1),1);

if options.show
    close all
end
w     = options.w;

r1 = 0.2*p;
r2 = 0.4*p;

i0 = p/2;
j0 = p/2;
Ra = zeros(p,p);
Rb = zeros(p,p);

for i=1:p
    for j=1:p
        r = sqrt((i-i0)^2+(j-j0)^2);
        if r<r2
            if r<r1
                Ra(i,j)=1;
            else
                Rb(i,j)=1;
            end
        end
    end
end

ra = Ra==1;
rb = Rb==1;


Z = zeros(5000,m);
id=0;
% XX = [];
for i=1:n
    Jas = zeros(p,p);    % mean of J
    Jhs = Jas;           % mean of Jh
    for j=1:m
        Ia  = Iw(p*(i-1)+1:p*i,p*(j-1)+1:p*j);
        Ih  = vl_harris(Ia,3);
        Jas = Jas + Ia;
        Jhs = Jhs + Ih;
    end
    Jas = Jas/m;
    Jhs = Jhs/m;
    Sa  = Jas-min(Jas(:));
    Sb  = Jhs;
    
    ga = mean(Sa(ra));
    gb = mean(Sa(rb));
    ha = mean(Sb(ra));
    % hb = mean(Sb(rb));
    
    x1 = abs(ga-gb)/gb;
    x2 = ha;
    % x3 = abs(ha-hb)/hb;
    %XX = [XX;x1 x2]
    if [x1 x2 1]*w<0
        id = id + 1;
        Z(id,1:ny) = Y(i,:);
    else
        nz(p*(i-1)+1:p*i) = 0;
        d(i) = 0;
    end
end
Wd = W(:,:,d==1);
Iz(nz==0,:)=[];
Z = Z(1:id,:);
z = sum(Z,1);
Z = Z(:,z>0);
% XX