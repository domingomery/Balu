% function [Fpq,Bpq] = Bsq_fundamentalSIFT(kp,p,q)
%
% Toolbox: Balu
%
%    Fundamental matrix between two views p and q of a sequence using SIFT
%    descriptors kp.
%
%  Example: case when keypoints must be extracted
%
%    f.path             = ''; % Balu directory as path or current directory
%    f.extension        = '.png';
%    f.prefix           = 'X';
%    f.digits           = 1;
%    f.gray             = 1;
%    f.subsample        = 1;
%    f.resize           = 0;
%    f.window           = [];
%    f.negative         = 0;
%    f.sequence         = 1:6;
%    f.imgmin           = 1;
%    f.imgmax           = 6;
%    options.show       = 1;
%    options.descriptor = 'harris+sift';
%    options.clean      = 0;
%    kp                 = Bsq_des(f,options);
%    F13                = Bsq_fundamentalSIFT(kp,1,3);
%    F14                = Bsq_fundamentalSIFT(kp,1,4);
%    figure(1); Bio_imgshow(f,1,[]); hold on
%    figure(2); Bio_imgshow(f,3,[]); hold on
%    figure(3); Bio_imgshow(f,4,[]); hold on
%    while(1)
%        figure(1);
%        disp('click a point in Figure 1...') % click
%        p = vl_click; m1 = [p(2) p(1) 1]';
%        plot(p(1),p(2),'g+')
%        figure(2)
%        Bmv_epiplot(F13,m1)
%        figure(3)
%        Bmv_epiplot(F14,m1)
%    end
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [Fpq,Bpq] = Bsq_fundamentalSIFT(kp,p,q)
img      = kp.img;
des      = kp.des;
fra      = kp.fra;
% Fundamental matrix for image p and q using SIFT matching and RANSAC
ip       = find(img==p);
iq       = find(img==q);
dp       = des(ip,:);
dq       = des(iq,:);
frap     = fra(ip,:);
fraq     = fra(iq,:);
[matches, scores] = vl_ubcmatch(dp', dq');
ii       = find(scores<17000);
xp       = frap(matches(1,ii),[2 1]);
xq       = fraq(matches(2,ii),[2 1]);
[Fpq, i] = Bmv_fundamentalRANSAC(xp',xq');

inl      = ii(i);
matchp   = matches(1,inl);
matchq   = matches(2,inl);
Bpq      = [ip(matchp) iq(matchq)];

