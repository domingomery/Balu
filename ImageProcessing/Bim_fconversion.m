% Bim_fconversion(fmt1,fmt2)
%
% Toolbox: Balu
%    Format conversion of images
%
%    It converts all images with format fmt1 into images with format fmt2.
%    fmt1 and fmt2 are strings.
%
% D.Mery, PUC-DCC, Jun 2010
% http://dmery.ing.puc.cl
%

function Bim_fconversion(fmt1,fmt2)


d = dir(['*.' fmt1]);
n = length(d);
m1 = length(fmt1);
m2 = length(fmt2);
if (n>0)
    for i=1:n
        s1 = d(i).name;
        s0 = s1;
        s0(end-m1+1:end) = [];
        s2 =[s0 fmt2];
        fprintf('%3d/%3d) converting image %s into %s...\n',i,n,s1,s2)
        I = imread(s1);
        imwrite(I,s2,fmt2);
    end
else
    error('There is no image with format %s.',fmt1);
end