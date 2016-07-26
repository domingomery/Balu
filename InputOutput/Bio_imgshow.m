function Bio_imgshow(f,i,map)
I = Bio_loadimg(f,i);
if size(I,3)==1;
    if max(I(:))>5
        I = I/255;
    end
end
if exist('map','var')
    imshow(Bio_loadimg(f,i),map)
else
    imshow(Bio_loadimg(f,i))
end