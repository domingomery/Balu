% function function [F,Bo] = Bsq_multifundamental(kp,options)
%
% Toolbox: Balu
%
%  (No calibrated) Multiple fundamental matrices of a sequence
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
%
%    op1.show           = 1;
%    op1.descriptor     = 'harris+sift';
%    op1.clean          = 0;
%    op2.img1           = 1;
%    op2.img2           = 4;
%    op2.m              = 3;
%    op2.show           = 1;
%    
%    kp                 = Bsq_des(f,op1);
%    Fo                 = Bsq_fmulti(kp,op2);
%    ii = and(Fo(:,1)==1,Fo(:,2)==3);
%    F13 = zeros(3,3); F13(:) = Fo(ii,3:11);
%    ii = and(Fo(:,1)==1,Fo(:,2)==4);
%    F14 = zeros(3,3); F14(:) = Fo(ii,3:11);
%
%    close all
%    figure(1); Bio_imgshow(f,1,[]); hold on
%    figure(2); Bio_imgshow(f,3,[]); hold on
%    figure(3); Bio_imgshow(f,4,[]); hold on
%    while(1)
%        figure(1);
%        disp('click a point in Figure 1...') % click
%        p = vl_click; m1 = [p(2) p(1) 1]';
%        plot(p(1),p(2),'g+')
%        figure(2)
%        Bmv_epiplot(F13,m1);
%        figure(3)
%        Bmv_epiplot(F14,m1);
%    end
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [F,Bo] = Bsq_multifundamental(kp,options)



img1 = options.img1;
img2 = options.img2;
show = options.show;
m    = options.m;
N    = max(kp.img);
F    = zeros(10000,11);
j    = 0;
Bo   = zeros(20000,2);
i    = 1;
if show
    disp('Computing Fundamental Matrices...')
end
for p = img1:img2
    for q = p+1:min([p+m N])
        j              = j+1;
        if show
            fprintf('F(%d,%d)...',p,q)
        end
        [Fpq,Bpq]      = Bsq_fundamentalSIFT(kp,p,q);
        F(j,:)         = [p q Fpq(:)'];
        in             = size(Bpq,1);
        Bo(i:i+in-1,:) = Bpq;
        i              = i+in;
    end
end
F = F(1:j,:);
Bo = Bo(1:i-1,:);
if show
    fprintf('\n');
end