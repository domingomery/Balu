% Bio_copyfiles(prefix1,prefix2)
%
% Toolbox: Balu
%    Copy files prefix1* to prefix2*.
%    
%    This program convert all files starting with prefix1 to new files with
%    prefix2.
%
% Example:
%    Bio_copyfiles('images','img') % converts all files 'images' to 'img'
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl
function Bio_copyfiles(prefix1,prefix2)

d1 = dir([prefix1 '*']);
t1 = length(prefix1);
n = length(d1);

if n>0
    fprintf('copying %d files...\n',n)
    if ispc
        fcp = '!copy ';
    else
        fcp = '!cp ';
    end
    
    for i=1:n
        f1 = d1(i).name;
        f2 = [prefix2 f1(t1+1:end)];
        eval([fcp f1 ' ' f2]);
    end
end
