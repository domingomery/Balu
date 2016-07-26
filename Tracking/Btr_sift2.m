% Bpq = Btr_sift2(kp,p,q,options)
%
% Toolbox: Balu
%
%    Matching points between views p and image q using SIFT keypoints.
%
%    kp keypoints structure according function Bsq_des (see help)
%    options.matching 1: matching is estimated directly using vl_ubcmatch
%    function.
%    options.matching 2: matching is estimated using vl_ubcmatch function
%    and homography correspondences.
%
%    Bpq is a Nx2 matrix with N matchings. Bpq(i,1) and Bpq(i,2) mean that
%    a matching between keypoint Bpq(i,1) and Bpq(i,2) was found.
%
%  Example:
%
%    f.path            = '/Volumes/domingomery/Mingo/Matlab/balu3/';
%                      %  ^^^          directory of Balu        ^^^
%    f.extension       = '.jpg';
%    f.prefix          = 'testimg';
%    f.digits          = 1;
%    f.gray            = 1;
%    f.subsample       = 1;
%    f.resize          = 0.5;
%    f.window          = [];
%    f.negative        = 0;
%    f.sequence        = 1:6;
%    f.imgmin          = 1;
%    f.imgmax          = 6;
%    options.matching  = 2;
%    options.show      = 1;
%    kp = Bsq_des(f,'harris+sift',options);
%    p = 5; q = 6;
%    B = Btr_sift2(kp,p,q,options);
%    figure
%    imshow(Bloadimg(f,p),[]); hold on; title('Image p')
%    for i=1:size(B,1)
%        plot(kp.fra(B(i,1),1),kp.fra(B(i,1),2),'*')
%        text(kp.fra(B(i,1),1),kp.fra(B(i,1),2),num2str(i))
%    end
%    figure
%    imshow(Bloadimg(f,q),[]); hold on; title('Image q')
%    for i=1:size(B,1)
%        plot(kp.fra(B(i,2),1),kp.fra(B(i,2),2),'*')
%        text(kp.fra(B(i,2),1),kp.fra(B(i,2),2),num2str(i))
%    end
%
%  See also Bsq_des, Btr_siftn.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function Bpq = Btr_sift2(kp,p,q,options)




method = options.matching;

img      = kp.img;
des      = kp.des;
fra      = kp.fra;
ip       = find(img==p);
iq       = find(img==q);
frap     = fra(ip,:);
fraq     = fra(iq,:);
dp       = des(ip,:);
dq       = des(iq,:);
[matches, scores] = vl_ubcmatch(dp', dq',1.1);


switch method

    case 1 % direct
        ii     = find(scores<100000); %era 30000
        xp     = frap(matches(1,ii),[2 1]);
        xq     = fraq(matches(2,ii),[2 1]);
        dx     = xp-xq;
        dx2    = sqrt(sum(dx.*dx,2));
        [i,j]  = sort(dx2);
        ii     = i<150;
        ni     = length(ii);
        j     = j(1:min([2000 ni]));
    case 2 % with homographies (adapted from Andrea Vedaldi)
        nm    = size(matches,2) ;
        X1    = frap(matches(1,:),1:2)' ; X1(3,:) = 1 ;
        X2    = fraq(matches(2,:),1:2)' ; X2(3,:) = 1 ;
        nt    = 100;
        score = zeros(nt,1);
        ok    = zeros(nt,nm);
        for t = 1:nt
            % estimate homography
            subset = vl_colsubset(1:nm, 4) ;
            A = [] ;
            for i = subset
                A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
            end
            [U,S,V] = svd(A) ;
            Ht = reshape(V(:,9),3,3) ;
            % score homography
            X2_ = Ht * X1 ;
            du  = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
            dv  = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
            ok(t,:) = (du.*du + dv.*dv) < 6*6 ;
            score(t) = sum(ok(t,:)) ;
        end
        [score, best] = max(score) ;
        j = ok(best,:)==1;

end
Bpq = [ip(matches(1,j)) iq(matches(2,j))];
