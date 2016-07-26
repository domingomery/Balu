% Bo = Btr_siftn(kp,options)
%
% Toolbox: Balu
%
%    Matching points in all two consecutive views of a sequence using SIFT 
%    keypoints.
%
%    kp keypoints structure according function Bsq_des (see help)
%    options.matching 1: matching is estimated using vl_ubcmatch
%    function only.
%    options.matching 2: matching is estimated using vl_ubcmatch function
%    and homography correspondences.
%    options.show display results.
%
%    Bo is a Nx2 matrix with N matchings. a=Bo(i,1) and b=Bo(i,2) mean that
%    a matching between keypoint a and b was found in views kp.img(a) and
%    kp.img(b).
%
%  Example:
%
%    f.path            = '/Volumes/domingomery/Mingo/Matlab/balu3/';
%                      %  ^^^          directory of Balu        ^^^
%    f.extension       = '.jpg';
%    f.prefix          = 'I';
%    f.digits          = 1;
%    f.gray            = 1;
%    f.subsample       = 1;
%    f.resize          = 0;
%    f.window          = [];
%    f.negative        = 0;
%    f.sequence        = [2 4 1 5 3 6];  % sorted sequence (see Btr_sort)
%    f.imgmin          = 1;
%    f.imgmax          = 6;
%    options.matching  = 2;
%    options.show      = 1;
%
%    kp = Bsq_des(f,'harris+sift',options);     % keypoints
%    Bo = Btr_siftn(kp,options);                % matching points
%    options.plottraj   = 1;                    % plot trajectories
%    options.plotimg    = 0;                    % display image numbers
%    options.plotpoints = 0;                    % display points numbers
%    options.plotsquare = 0;                    % display points numbers
%    Btr_plot(kp,Bo,f,options)                  % plot matched tracks
%
%  See also Bsq_des, Btr_sift2.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function Bo = Btr_siftn(kp,options)

show = options.show;

j  = 0;
Bo = zeros(20000,2);
i  = 1;
p1 = min(kp.img);
p2 = max(kp.img);

for p = p1:p2-1
    q = p+1;
    j = j+1;
    if show
        fprintf('Btr_siftN  : Matching (%d,%d)...\n',p,q)
    end
    Bpq = Btr_sift2(kp,p,q,options);
    in  = size(Bpq,1);
    Bo(i:i+in-1,:) = Bpq;
    i   = i+in;
end
Bo = Bo(1:i-1,:);
n = size(Bo,1);
if show
    fprintf('Btr_siftN  : %4d matchings in the whole sequence.\n',n)
end