% [Z,f] = Btr_analysis(kp,Y,P,options)
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
%    to function Bio_loadimg (see help).
%
%
%  Example:
%      See example in Btr_demo.
%
%  See also Bsq_des, Btr_join, Btr_merge, Btr_2, Btr_3, Btr_demo.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [W,Iw] = Btr_windows(kp,Y,P,files,options)
show  = options.show;
if show
    close all
end

if isfield(options,'rotate');
    rot = options.rotate;
else
    rot = 0;
end


n   = size(Y,1);
m   = max(kp.img);
fra = kp.fra;
img = kp.img;
s   = 1.7;
p   = 128;
p2  = round(p/2);
W   = zeros(m,4,n);

I   = Bio_loadimg(files,files.imgmin);
[Ni,Mi] = size(I);
Ij = zeros(Ni,Mi,m);
Ij(:,:,1) = I;
for j=2:m
    Ij(:,:,j) = Bio_loadimg(files,files.imgmin+j-1);
end
Iw = zeros(p*n,p*m);

for i=1:n
    Yi = Y(i,:);
    Yn = Yi(Yi>0);
    myi = length(Yn);
    mm = zeros(3,myi);
    Ps = zeros(3*myi,4);
    for k=1:myi
        mm(:,k) = [double(fra(Yn(k),[2 1])) 1]';
        im = img(Yn(k));
        Ps(3*k-2:3*k,:) = P(im*3-2:im*3,:);
    end
    M = Bmv_reco3dn(mm,Ps);
    r  = mean(fra(Yn,3));
    b  = mean(fra(Yn,4));
    for j=1:m
        [t,tj] = max(j==img(Yn));
        if t % region was detected we use the original coordinates
            a  = Yn(tj);
            x  = fra(a,1);
            y  = fra(a,2);
            r  = fra(a,3);
            b  = fra(a,4);
        else % region was not detected we use the reprojection
            mj = P(j*3-2:j*3,:)*M; mj = mj/mj(3);
            x = mj(2);
            y = mj(1);
            % r is the last radius or the average
        end
        
        if rot
            i1 = max([1 round(y-2*s*r)]); i2 = min([round(y+2*s*r) Ni]);
            j1 = max([1 round(x-2*s*r)]); j2 = min([round(x+2*s*r) Mi]);
            i10 = max([1 round(y-s*r)]); i20 = min([round(y+s*r) Ni]);
            j10 = max([1 round(x-s*r)]); j20 = min([round(x+s*r) Mi]);
            Io = imresize(Ij(i1:i2,j1:j2,j),[2*p 2*p]);
            Io = imrotate(Io,90+b*180/pi,'bilinear','crop');
            Io = Io(p2:3*p2-1,p2:3*p2-1);
        else
            %i10 = max([1 round(y-s*r)]); i20 = min([round(y+s*r) Ni]);
            %j10 = max([1 round(x-s*r)]); j20 = min([round(x+s*r) Mi]);
            i10 = max([1 min([round(y-s*r) Ni])]); i20 = min([round(y+s*r) Ni]);
            j10 = max([1 min([round(x-s*r) Mi])]); j20 = min([round(x+s*r) Mi]);
            % W(j,:,i) = [i10 i20 j10 j20];
            % [i10 i20 j10 j20 size(Ij)]
            Io = imresize(Ij(i10:i20,j10:j20,j),[p p]);
        end
        W(j,:,i) = [i10 i20 j10 j20];
        Iw(p*(i-1)+1:p*i,p*(j-1)+1:p*j) = Io;
    end
end

if show
    figure(1)
    for j=1:m
        if j==1
            clf
            imshow(Ij(:,:,j),[])
            hold on
        end
        for i=1:n
            i1=W(j,1,i);
            i2=W(j,2,i);
            j1=W(j,3,i);
            j2=W(j,4,i);
            plot([j1 j2 j2 j1 j1],[i1 i1 i2 i2 i1],'r');
        end;
    end
    figure(2)
    imshow(Iw,[])
end

