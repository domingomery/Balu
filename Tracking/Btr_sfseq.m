% [P,f,kp,H1] = Btr_sfseq(f,options)
%
% Toolbox: Balu
%
%    Structure from an image sequence
%
%    The structure estimation is obtained by computing a geometric model of
%    the multi-views from the object itself. The geometric model is
%    estimated by a bundle adjustment algorithm on stable SIFT keypoints
%    across multi-views that are not necessary sorted.
%
%    f is the structure that defines the sequence
%
%    f.path          : directory where are the files
%    f.extension     : extension (eg: 'jpg')
%    f.prefix        : prefix (eg: 'DSC_')
%    f.digits        : number of digits (eg:4)
%    f.gray          : 1 means rgb to gray conversion
%    f.subsample     : subsampling rate (eg:1 means no subsample)
%    f.resize        : parameter of imresize, 0 means no imresize
%    f.window        : image window
%    f.negative      : negative window
%    f.sequence      : if seq = [3 2 1], the image for i=1 will be No. 3
%    f.imgmin        : minimal number of the images (not used by Bloadimg)
%    f.imgmax        : maximal number of the images (not used by Bloadimg)
%
%    options.show         1: display results
%    options.descriptor   keypoints' descriptors
%    options.sfm_method   SfM > 1: affine, 2: projective (see Btr_sfm)
%    options.sfm_samples  all samples for RANSAC SfM (see Btr_sfm)
%    options.sfm_iter     iterations for RANSAC SfM (see Btr_sfm)
%    options.sfm_cover    cover of the SfM solution (see Btr_sfm)
%    options.matching     1: direct, 2: homomography (see Btr_sift2)
%    options.indexonly    sort and change image sequence (see Bsq_sort)
%    options.sort         1: sort the image sequence
%
% Example:
% f.path          = '/Volumes/domingomery/Mingo/Matlab/balu3/';
% f.prefix        =  'X';   f.extension   =  '.png';
% f.digits        = 1;
% f.gray          = 1;
% f.subsample     = 1;
% f.resize        = 1;
% f.imgmin        = 1;
% f.imgmax        = 6;
% f.window        = [];
% f.negative      = 0;
% f.sequence      = 1:f.imgmax;
%
% show = 0;                       % use show = 1 to display results
% op1.show        = show;         % 1: display results
% op1.descriptor  = 'harris';     % keypoints' descriptors
% op1.sfm_method  = 1;            % SfM > 1: affine, 2: projective (see Btr_sfm)
% op1.sfm_samples = 0;            % all samples for RANSAC SfM (see Btr_sfm)
% op1.sfm_iter    = 1;            % iterations for RANSAC SfM (see Btr_sfm)
% op1.sfm_cover   = 0.5;          % cover of the SfM solution (see Btr_sfm)
% op1.matching    = 1;            % 1: direct, 2: homomography (see Btr_sift2)
% op1.sort        = 1;            % sort image sequence using Bsq_sort
% op1.indexonly   = 0;            % sort and change image sequence (see Bsq_sort)
%
% [P,f,kp1] = Btr_sfseq(f,op1);
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl



function [P,f,kp1,H1] = Btr_sfseq(f,options)


nimg              = f.imgmax-f.imgmin+1;  % number of images in the sequence
options.join_iter = nimg-2;               % number of join iterations (see Btr_join)


kp1     = Bsq_des(f,options);             % Keypoints

if options.sort
    [kp1,f] = Bsq_sort(kp1,f,options);    % Sorted images
end

Bo          = Btr_siftn(kp1,options);     % Matching Points

if nimg>2
    Ho      = Btr_join(Bo,[],options);    % structure tracks
else
    Ho      = Bo;
end

ok = 0;
while not(ok)
    [P,H1]      = Btr_sfm(kp1,Ho,options);    % Bundle Adjsutment
    if isempty(H1)
        options.sfm_samples = options.sfm_samples+4;
        options.sfm_iter    = options.sfm_iter+50;
    else
        ok=1;
    end
end

