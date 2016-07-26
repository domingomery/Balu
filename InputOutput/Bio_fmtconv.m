% Bio_fmtconv(fmt1,fmt2)
%
% Toolbox: Balu
%    Image format conversion from format fmt1 to fmt2.
%    
%    This program convert all fmt1 images of current directory to fmt2
%    images. fmt1 and fmt2 are strings.
%
% Example:
%    Bio_fmtconv('jpg','png') % converts all jpg images into png images
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl
function Bio_fmtconv(fmt1,fmt2)

f1 = dir(['*.' fmt1]);
t1 = length(fmt1);
n = length(f1);
if n>0
    for i=1:n
        fi1 = f1(i).name;
        I = imread(fi1);
        fi2 = [fi1(1:end-t1) fmt2];
        imwrite(I,fi2,fmt2);
    end
end
