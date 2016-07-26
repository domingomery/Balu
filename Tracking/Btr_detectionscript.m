% [Tf,kp1,kp2,f,P] = Btr_detection(f,op1,op2,op3)
%
% Toolbox: Balu
%
%    Detection by tracking.
%    Fine Detection in Rigid Objects using Multi-views
%
%    This method consists of two steps: Structure Estimation, to obtain a
%    geometric model of the multi-views from the object itself, and
%    Details Detection, to detect the relevant details of the object.
%
%    The geometric model is estimated by a bundle adjustment algorithm on
%    stable SIFT keypoints across multi-views that are not sorted. The
%    details detection is performed by a extremely simple segmentation
%    algorithm followed by a tracking algorithm based on geometric and
%    appearance constraints. It is not required that the relevant details
%    have to be segmented in all views. Additionally, it is allowed to
%    obtain false detections in this step. The tracking is used to
%    eliminate the false detections without discriminating the relevant
%    details.
%
% Example 1: Uncalibrated approach (structure must be estimated)
% *** see the code of Btr_demo
%
% Example 2: Calibrated approach, projection matrices are already estimated
% % 1. Structure
% % Store in matrix P the projection matrices as follows P(k*3-2:k*3,:) is
% % the 3x4 porjection of view k, i.e., P = [P1;P2;P3;P4;P5;P6];
% % *** if you want to estimate P using a SfM approache see Btr_sfseq ***
% % *** using P from last example the following steps can be used after the
% % *** the option definitions of Example 1:
% % *** op1.sort = 0;
% % *** [P,f,kp1] = Btr_sfseq(f,op1);
% % *** after estimation of P, we can use:
%
% % 2. Details Detection & Display final results
% [Tf,kp1,kp2,f,P] = Btr_detection(f,P,op2,op3);
%
% See also Btr_demo, Btr_gui.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


% function [Tf,kp1,kp2,f,P,Iz,Wd,inx] = Btr_detection(f,op1,op2,op3)

% 0. Loading images
% =================
f    = Bsq_load(f);
nimg = f.imgmax-f.imgmin+1;


% 1. Structure Estimation
% =======================
if sum(class(op1)=='struct')==6
    [P,f,kp1,H1] = Btr_sfseq(f,op1);
else
    P = op1;
    kp1 = [];
    H1 = [];
end


% 2. Details Detection
% =====================
if isfield(op2,'kp')
    kp2 = op2.kp;
else
    kp2 = Bsq_des(f,op2);               % Segmentation+Description
end
F       = Bsq_fundamental(P);           % Fundamental Matrices
T2      = Btr_2(kp2,F,op2);             % Tracking in 2 views

if nimg>2
    T       = Bsq_trifocal(P);          % Trifocal Tensors
    T3 = Btr_3(kp2,T2,T,op2);           % Tracking in 3 views
    if nimg>3
        op2.join_iter = 1;
        Tr = Btr_join(T3,[],op2);       % Matching in 4 views
    else
        Tr = T3;
    end
else
    Tr = T2;
end
Tm        = Btr_merge(kp2,Tr);          % Merging tracks
[W,Iw]    = Btr_windows(kp2,Tm,P,f,op2);% Multi-view windows (MVW)
[Tf,Iz,Wd]= Btr_classify(W,Iw,Tm,op2);  % Classification of MVW

% 3. Display results
% ====================
if exist('op3','var')
    if ~isfield(op3,'show')
        op3.show = 1;
    end
    if op3.show
        figure(2)
        op3.img1        = 2;
        Btr_plot(kp2,Tf,f,op3)
        pause(1)
        figure(1)
        op3.img1        = 1;
        Btr_plot(kp2,Tf,f,op3)
        figure(3)
        imshow(Iz,[])
        fprintf('%3d regions detected.\n',size(Tf,1))
    end
end

inx.H1 = H1;
inx.T2 = T2;
if nimg>2
    inx.T3 = T3;
else
    inx.T3 = T2;
end
inx.Tr = Tr;
inx.Tm = Tm;
inx.Tf = Tf;
