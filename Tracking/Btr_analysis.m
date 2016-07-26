% [Z,f] = Btr_analysis(kp,Y,P,files,options)
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

function [Z,f] = Btr_analysis(kp,Y,P,files,options)

if options.show
    close all
end
[n,m] = size(Y);
nimg    = files.imgmax-files.imgmin+1;
fra   = kp.fra;
img   = kp.img;
show  = options.show;
w     = options.w;
sc    = 0.25;

projective = 0;

s = 1.7;
p = 128;

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
f = zeros(5000,3);
id=0;
JJ = [];
for i=1:n
    Yi = Y(i,:);
    Yn = Yi(Yi>0);
    myi = length(Yn);
    m = zeros(3,myi);
    Ps = zeros(3*myi,4);
    for k=1:myi
        m(:,k) = [double(fra(Yn(k),[2 1])) 1]';
        im = img(Yn(k));
        Ps(3*k-2:3*k,:) = P(im*3-2:im*3,:);
    end
    [M,err] = Bmv_reco3dn(m,Ps);
    if nimg>=4
        me = conv2(err,ones(1,4)/4,'valid');
    else
        me = err;
    end
    if min(me)<7
        % im  = ones(nimg,1);
        Ja  = zeros(p,nimg*p); % detected windows
        Jh  = Ja;            % harris of J
        Jas = zeros(p,p);    % mean of J
        Jhs = Jas;           % mean of Jh
        r  = mean(fra(Yn,3));
        II = zeros(p,p,nimg);
        Is = [];
        ns = 0;
        xx = zeros(nimg,1);
        yy = xx;
        rr = xx;
        xj = xx;
        for j=1:nimg
            I = Bio_loadimg(files,files.imgmin+j-1);
            [Ni,Mi] = size(I);
            [t,tj] = max(j==img(Yn));
            if t % region was detected we use the original coordinates
                a  = Yn(tj);
                x  = fra(a,1);
                y  = fra(a,2);
                r  = fra(a,3);
            else % region was not detected we use the reprojection
                mj = P(j*3-2:j*3,:)*M; mj = mj/mj(3);
                x = mj(2);
                y = mj(1);
                % r is the last radius or the average
            end
            if show
                figure(1)
                clf
                imshow(I,[])
                vl_plotframe([x y r]','g','LineWidth',1);
                Is = [Is imresize(I,sc)]; %#ok<AGROW>
                
                enterpause(0)
            end
            xj(j) = x+(j-1)*Mi;
            xx(j) = x;
            yy(j) = y;
            rr(j) = r;
            i1 = round(y-s*r); i2 = round(y+s*r);
            j1 = round(x-s*r); j2 = round(x+s*r);
            if (i1>0) && (i2<=Ni) && (j1>0) && (j2<=Mi)
                Iw                  = I(i1:i2,j1:j2);
                Ia                  = imresize(Iw,[p p]);
                Ja(:,(j-1)*p+1:j*p) = Ia;
                Ih                  = imresize(vl_harris(Iw,3),[p p]);
                Jh(:,(j-1)*p+1:j*p) = Ih;
                II(:,:,j)           = Ia;
                Jas                 = Jas + Ia;
                Jhs                 = Jhs + Ih;
                ns = ns + 1;
            end
        end
        if projective
            ni = 0;
            HH = zeros(3,3,nimg-1);
            kh = 0;
            for j=2:nimg
                Ij = II(:,:,j);
                Ik = II(:,:,j-1);
                [Iks,Hjk,err]=Bmv_homography(Ij,Ik,0);
                if err==1
                    HH(:,:,j-1) = Hjk;
                    kh = kh+1;
                end
            end
            if kh==nimg-1
                JJ = zeros(p,p);
                for j=1:nimg-1
                    for k=j:nimg-1
                        H  = HH(:,:,k);
                        Ik = II(:,:,j);
                        Iks = Bmv_projective2D(Ik,H,size(Ik),0);
                        II(:,:,j) = Iks;
                    end
                end
                JJJ = [];
                for j=1:nimg
                    JJ = JJ+II(:,:,j);
                    JJJ = [JJJ II(:,:,j)]; %#ok<AGROW>
                end
                Js = JJ/nimg;
                if show
                    figure(3)
                    imshow([[Ja Js/ni];[JJJ JJ]],[])
                    enterpause
                end
            end
            x1 = 0;
            x2 = 0;
            x3 = 0;
            % J = JJJ;
        else
            Jas = Jas/ns;
            Jhs = Jhs/ns;
            Sa  = Jas-min(Jas(:));
            Sb  = Jhs;
            %ta  = mean(mean(Sa(P1:P2,P1:P2)))/mean(mean(Sa(P3:P4,P3:P4)));
            %tb  = mean(mean(Sb(P1:P2,P1:P2)));
            
            ga = mean(Sa(ra));
            gb = mean(Sa(rb));
            ha = mean(Sb(ra));
            hb = mean(Sb(rb));
            
            x1 = abs(ga-gb)/gb;
            x2 = ha;
            x3 = abs(ha-hb)/hb;
            % [x1 x2 x3]
        end
        if [x1 x2 1]*w<0
            %            if x1*w(1)
            
            
            id = id + 1;
            Z(id,1:myi) = Yn;
            f(id,:) = [x1 x2 x3];
            if show
                Ji = [Bimglin(Ja) Bimglin(Jas) Bimglin(Jhs)];
                JJ = [JJ;Ji]; %#ok<AGROW>
                figure(3)
                vl_plotframe([xx yy rr]','r','LineWidth',1);
                plot(xx,yy)
                figure(1)
                clf
                imshow(Is,[])
                vl_plotframe(sc*[xj yy 3*rr]','r','LineWidth',1);
                rg = get(1,'OuterPosition'); rgn1 = [1 500 rg(3) rg(4)];set(1,'OuterPosition',rgn1)
                figure(2)
                imshow(Ji,[])
                rg = get(2,'OuterPosition'); rgn1 = [1 200 rg(3) rg(4)];set(2,'OuterPosition',rgn1)
                % enterpause
                
                
            end
        end
    end
end
Z = Z(1:id,:);
z = sum(Z,1);
Z = Z(:,z>0);
f = f(1:id,:);
if show
    figure(4)
    imshow(JJ,[])
    figure(3)
end


