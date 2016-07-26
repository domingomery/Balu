% I = Bloadimg(f,i)
%
% Toolbox: Balu
%
% Load image i of set image defined by structure f
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
%    f.show          : 1 means display image (deafult = 0)
%
%    Example:
%    f.path          = ''; % Balu directory as path or current directory
%    f.extension     = '.jpg';
%    f.prefix        = 'testimg';
%    f.digits        = 1;
%    f.gray          = 1;
%    f.subsample     = 1;
%    f.resize        = 0;
%    f.window        = [];
%    f.negative      = 0;
%    f.sequence      = 1:6;
%    f.imgmin        = 1;
%    f.imgmax        = 6;
%    I = Bloadimg(f,3);
%    imshow(I,[])
%
%    See Bseq_show.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [I,st] = Bio_loadimg(f,i)


ipath = f.path;
iext  = f.extension;
ipre  = f.prefix;


if isfield(f,'resize')
    re     = f.resize;
else
    re = 0;
end

if isfield(f,'subsample')
    t     = f.subsample;
else
    t = 1;
end

if isfield(f,'show')
    show    = f.show;
else
    show = 0;
end

if isfield(f,'window')
    w    = f.window;
else
    w = [];
end


if isfield(f,'negative')
    neg    = f.negative;
else
    neg = 0;
end

if isfield(f,'gray')
    gr    = f.gray;
else
    gr = 0;
end

if isfield(f,'sequence')
    seq    = f.sequence;
else
    seq = f.imgmin:f.imgmax;
end


if ~isempty(ipath)
    if ipath(end)~='/'
        ipath = [ipath '/'];
    end
end

%try
%    I = f.images(:,:,f.sequence(i)-f.imgmin+1);
%
%catch exception
if isfield(f,'images')
    I = f.images(:,:,f.sequence(i)-f.imgmin+1);
    % I = f.images(:,:,i-f.imgmin+1);
    st = 'from memory';
    
else
    if show
    fprintf('Loading image %d...',i);
    end
    if ipre == '*'
        if iext(1)=='.'
            iext = iext(2:end);
        end
        dfiles = [dir([ipath '*.' lower(iext)]); dir([ipath '*.' upper(iext)])];
        %tf = cat(1,dfiles.name);
        %[tj,ti] = sort(tf);
        ti = 1:length(dfiles);
        if (i>0)&&(i<=length(dfiles))
            st = [ipath dfiles(ti(i)).name];
            if ~exist(st,'file')
                error('file %s does not exist.',ipath)
            end
            if compare(upper(iext),'.MAT')==0
                st(end-2:end) = 'mat';
                s = ['load ' st];
                eval(s);
            else
                I = imread(st);
                if show
                    disp(st)
                end
            end
        else
            st = -1;
            if i==0
                I = length(tf);
            else
                I = -1;
            end
        end
    else
%        if iext(1)~='.'
%            iext = ['.' iext]
%        end
        d     = f.digits;
        sti = sprintf('0000000%d',seq(i));
        sti = sti(end-d+1:end);
        st = sprintf('%s%s%s%s',ipath,ipre,sti,iext);
        if (exist(st,'file'))
            if compare(upper(iext),'.MAT')==0
                st(end-2:end) = 'mat';
                s = ['load ' st];
                eval(s);
            else
                I = imread(st);
                if show
                    disp(st)
                end
            end
        else
            I = -1;
            fprintf('%s does not exist.\n',st)
        end
    end
    
    if length(I(:))>1
        if ~isempty(w)
            I = I(w(1):w(2),w(3):w(4),:);
        end
        if size(I,3)==3
            if gr
                I = rgb2gray(I(1:t:end,1:t:end,:));
            else
                I = I(1:t:end,1:t:end,:);
            end
        else
            I = I(1:t:end,1:t:end,:);
        end
        if neg
            I = 256.0-double(I);
        end
    end
    I = double(I);
    if re>0
        I = imresize(I,re);
    end
    if show
        imshow(I/256,[])
    end
end