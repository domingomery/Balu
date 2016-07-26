% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function d = Bio_labelimage(f)


d = zeros(f.imgmax-f.imgmin+1,1);

for i=f.imgmin:f.imgmax
    [I,st] = Bio_loadimg(f,i);
    imshow(I(:,:,1),[])
    fprintf('\n--- processing image %s...\n',st);
    d(i-f.imgmin+1,1) = input('Label for this image? ');
end