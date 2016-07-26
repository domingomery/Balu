% [rk,j] = Bsq_vgoogle(f,i,ilu,v,show)
%
% Toolbox: Balu
%
%    Search sequence images similar to image i.
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
%    f.imgmin        : minimal number of the images (not used by Bio_loadimg)
%    f.imgmax        : maximal number of the images (not used by Bio_loadimg)
%
%    ilu is a look up table of the real number of the images (see example),
%    if ilu is empty ilu will be f.imgmin:f.imgmax.
%    v is the document representation using "term frequency-inverse
%    document frequency"
%    show display results
%
%    Reference:
%      Sivic & Zisserman: Efficient visual search for videos cast as text
%      retrieval. 31(4):591-606. PAMI 2009.
%
%    Example:
%    f.path             = ''; % Balu directory as path or current directory
%    f.extension        = '.jpg';
%    f.prefix           = 'testimg';
%    f.digits           = 1;
%    f.gray             = 1;
%    f.subsample        = 1;
%    f.resize           = [256 256];
%    f.window           = [];
%    f.negative         = 0;
%    f.sequence         = 1:6;
%    f.imgmin           = 1;
%    f.imgmax           = 6;
%    options.show       = 1;
%    options.descriptor = 'sift';
%    options.clean      = 0;
%    kp = Bsq_des(f,options);
%    v  = Bsq_vocabulary(kp,100,options);
%    [rk,j] = Bsq_vgoogle(f,5,[],v,1);  % similar images to image 5
%
% See also Bsq_sort, Bsq_vocabulary.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [rk,j] = Bsq_vgoogle(f,i,ilu,v,show)


if isempty(ilu)
    ilu = f.imgmin:f.imgmax;
end

vq = v(ilu==i,:)';
[rk,j] = Bfa_vecsimilarity(vq,v);

if ~exist('show','var')
    show = 0;
end

if show
    disp('Finding similar images:')
    close all
    figure(1)
    imshow(Bio_loadimg(f,i),[])
    title(sprintf('query image %d',i))
    m = 4;
    for k=m+1:-1:2
        figure(k)
        imshow(Bio_loadimg(f,ilu(j(k))),[]);
        s = sprintf('Figure %d: image %d (Similarity with image %d: %f)',k,ilu(j(k)),i,rk(k));
        title(s)
        disp(s)
    end
    figure(1)
    disp('Figure 1: Query image')
end

