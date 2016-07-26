% f_new = Bsq_load(f)
%
% Toolbox: Balu
%
%    Load an image sequence.
%    
%    f is a file structure (see Bio_loadimg for details)
%  
%    f_new includes a fild called f_new.images with an array NxMxm for a 
%    sequence with NxM images
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
%    f = Bsq_load(f);
%    Ij = f.images;
%    imshow(Ij(:,:,3),[])
%
% See also Bload_img.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


function f_new = Bsq_load(f)

nimg = f.imgmax-f.imgmin+1;

I = Bio_loadimg(f,f.imgmin);
[N,M] = size(I);
Ij = zeros(N,M,nimg);
Ij(:,:,1) = I;
for j=f.imgmin+1:f.imgmax
    Ij(:,:,j-f.imgmin+1) = Bio_loadimg(f,j);
end
f_new = f;
f_new.images = Ij;
