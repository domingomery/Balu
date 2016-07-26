% Btr_demo
%
% Toolbox: Balu
%
%    Fine Detection in Rigid Objects using Multi-views
%    
%    In this demo we detect pen tips in X-ray images of a pencil case using
%    multiviews.
%
%    This demo consists of two steps: Structure Estimation, to obtain a
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
%    In this demo we use Btr_detection that calls the following functions:
%
%    Bsq_des             + Description of a sequence 
%    Bsq_load            + Load an image sequence
%    Bsq_show            + Display of a sequence
%    Bsq_sort            + Sort of a sequence 
%    Bsq_vocabulary      + Visual vocabulary of a sequence 
%    Bsq_fundamental     + Fundamental matrices of a sequence
%    Bsq_trifocal        + Trifocal tensors of a sequence
%    Btr_2               + Matching in 2 views with epipolar cosntraint
%    Btr_3               + Matching in 3 views with trifocal cosntraint
%    Btr_join            + Join of matching points
%    Btr_merge           + Merge tracks with common matching points
%    Btr_plot            + Plot of tracks
%    Btr_sfm             + Structure from Motion
%    Btr_sfseq           + Structure from an image sequence
%    Btr_sift2           + Matching keypoints in 2 views of a sequence
%    Btr_siftn           + Matching keypoints in n views of a sequence
%    Btr_windows         + Multi-view windows (MVW)
%    Btr_classify        + Track classify using Multi-view windows
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


close all
clear all
clc
warning off %#ok<WNOFF>


disp('Demo Fine Detection in Rigid Objects using Multi-views')
disp(' ');
disp('Detection of pen tips in X-ray images of a pencil case')
disp(' ');
% show = input('display results? (1=yes, 0 =no) ');
show = 0;
% 0. Initial data
% ===============
f.path          = '';  % current directory or a path directory
f.prefix        =  'X';   f.extension   =  '.png';
f.digits        = 1;
f.gray          = 1;
f.imgmin        = 1;
f.imgmax        = 6;
f.sequence      = 1:f.imgmax;

% 1. Structure Estimation
% =======================
op1.show        = show;           % 1: display results
op1.descriptor  = 'harris';       % keypoints' descriptors
op1.sfm_method  = 1;              % SfM > 1: affine, 2: projective (see Btr_sfm)
op1.sfm_samples = 0;              % all samples for RANSAC SfM (see Btr_sfm)
op1.sfm_iter    = 1;              % iterations for RANSAC SfM (see Btr_sfm)
op1.sfm_cover   = 0.5;            % cover of the SfM solution (see Btr_sfm)
op1.matching    = 1;              % 1: direct, 2: homomography (see Btr_sift2)
op1.sort        = 1;              % sort image sequence using Bsq_sort
op1.indexonly   = 0;              % sort and change image sequence (see Bsq_sort)
op1.clean       = 0;


% 2. Details Detection
% ====================
op2.show        = show;           % 1: display results
op2.descriptor  = 'blops';        % Segmentation
op2.param       = [13 20 0 50 ...
    15 200 0.5 1.2];              % parameters for segblops
                                  % Matching in two views
op2.dxi2max     = 1000;           % Xi^2 distance between two descriptors
op2.dnmax       = 60;             % normalized distance
op2.mviews      = 3;              % consecutive views
op2.depimax     = 15;             % epipolar distance
op2.dtrimax     = 25;             % trifocal reproyection
op2.w           = [-8.75 -1 3]';  % parameters of analysis
op2.clean       = 0;


% 3. Display results
% ==================
op3.plottraj    = 1;              % plot trajectories
op3.plotimg     = 0;              % display image numbers
op3.plotpoints  = 0;              % display points numbers
op3.plotsquare  = 1;              % display squares

tic
[Tf,kp1,kp2,f,P,Iz,Wd,inx] = Btr_detection(f,op1,op2,op3);
toc
