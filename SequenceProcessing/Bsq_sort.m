% s = Bsq_sort(kp,f,options)
% [kp_new,files_new] = Bsq_sort(kp,files,options)
%
% Toolbox: Balu
%
%    Sort an image sequence.
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
%    kp are the keypoints of the sequence (extracted with Bsq_des or any
%    similar procedure). If kp is empty the kypoints will be extracted with
%    Bsq_des using 'harris+sift' method (see help Bsq_des).
%
%    options.show display results
%    options.indexonly 1: returns only s
%    options.indexonly 0: returns kp_new and files_new, they are the
%    corresponding keypoints and files sequence for the sorted images
%
%    s is the sorted sequence.
%
%  Example 1: case when keypoints must be extracted
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
%    options.indexonly  = 0;
%    options.descriptor = 'harris+sift';
%    options.clean      = 0;
%    [kp,f2] = Bsq_sort([],f,options);
%    Bsq_show(f); title('unsorted sequence')
%    figure
%    Bsq_show(f2); title('sorted sequence')
%
%
%  Example 2: case when keypoints are already extracted
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
%    options.indexonly  = 0;
%    options.descriptor = 'harris+sift';
%    kp = Bsq_des(f,options);
%    figure
%    Bsq_show(f); title('unsorted sequence')
%    [kp2,f2] = Bsq_sort(kp,f,options);
%    figure
%    Bsq_show(f2); title('sorted sequence')
%
%  Example 3: only the new indices are riquired
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
%    options.show       = 0;
%    options.indexonly  = 1;
%    options.descriptor = 'harris+sift';
%    kp = Bsq_des(f,options);
%    s = Bsq_sort([],f,options)
%
%  See also Bsq_vgoogle, Bsq_des.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [kp_new,files_new] = Bsq_sort(kp,files,options)

f = files;

show = options.show;

if isempty(kp)
    kp = Bsq_des(f,options);
end
v = Bsq_vocabulary(kp,200,options);

n = max(kp.img);

H = zeros(n,n);

for i=1:n
    q = kp.ilu==(i+f.imgmin-1);
    vq = v(q,:)';
    [rk,j] = Bfa_vecsimilarity(vq,v);
    H(i,j) = rk;
end

spmax = 0;

for j=1:n
    i0 = j;
    sj    = zeros(n,1);
    sj(1) = i0;
    w    = ones(n,1);
    w(i0) = 0;
    sp = 0;
    for i=2:n
        h     = w.*H(:,i0);
        [p,q] = max(h);
        i0    = q;
        sj(i)  = q;
        w(q)  = 0;
        sp = sp+p;
    end
    if sp>spmax
        spmax = sp;
        s    = sj;
    end

end


if show
    figure
    Bsq_show(f);title('unsorted sequence')
    f.sequence(f.imgmin:f.imgmax) = s+f.imgmin-1;
    figure
    Bsq_show(f);title('sorted sequence')
end


if options.indexonly
    kp_new = s;
    files_new = (sum(H)-1)/(n-1); % confidence, small values are not similar enough
else

    kp_new = kp;
    n = max(kp.img);
    for i=1:n
        ii = kp.img==s(i);
        kp_new.img(ii) = i;
    end

    files_new = f;

    files_new.sequence(files.imgmin:files.imgmax) = s+files.imgmin-1;
end



