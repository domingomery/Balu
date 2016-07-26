% Iseq = Bsq_show(f,n,map,k)
%
% Toolbox: Balu
%
% Display an image sequence defined by structure f.
%
%    n is the number of images per row (default is the number of images of
%    the sequence).
%
%    map is the map of the image, if not given will be used "[]".
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
%    Iseq is the output image with all images of the sequence.
%
%    Example:
%    f.path          = ''; % Balu directory as path or current directory
%    f.extension     = '.jpg';
%    f.prefix        = 'testimg';
%    f.digits        = 1;
%    f.gray          = 1;
%    f.subsample     = 1;
%    f.resize        = [100 100];
%    f.window        = [];
%    f.negative      = 0;
%    f.sequence      = 1:6;
%    f.imgmin        = 1;
%    f.imgmax        = 6;
%    Bsq_show(f,3);
%
%    See Bio_loadimg.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function Iseq = Bsq_show(f,n,map,k)

if (~exist('n','var'))
    n = f.imgmax-f.imgmin+1;
end

t = 0;
II = [];
Iseq = [];
for i=f.imgmin:f.imgmax
    s = f.sequence(i);
    s = i;
    Ii = Bio_loadimg(f,s);
    Iseq = [Iseq Ii];
    t = t + 1;
    if (t==n)
        II = [II; Iseq];
        t = 0;
        Iseq = [];
    end
end
if t>0
    Iseq = [Iseq zeros(size(Ii,1),size(Ii,2)*(n-t))];
    II = [II; Iseq];
end

Iseq = II;
if exist('map','var')
    if exist('k','var')
        if k>0
            imshow(k*II,map)
        end
    else
        imshow(II,map)
    end
else
    imshow(II,[])
end