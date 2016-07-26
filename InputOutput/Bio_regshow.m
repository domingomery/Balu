function Bio_regshow(I,R,c)

[N,M] = size(I);
Ic = zeros(N,M,3);
I = uint8(I);
X = I.*uint8(not(R));
Y = I;
Y(R==1)=255;
switch c
    case 1
        Ic(:,:,1)=X;
        Ic(:,:,2)=I;
        Ic(:,:,3)=I;
    case 2
        Ic(:,:,1)=I;
        Ic(:,:,2)=X;
        Ic(:,:,3)=I;
    case 3
        Ic(:,:,1)=I;
        Ic(:,:,2)=I;
        Ic(:,:,3)=X;
    case 4
        Ic(:,:,1)=Y;
        Ic(:,:,2)=I;
        Ic(:,:,3)=I;
    case 5
        Ic(:,:,1)=I;
        Ic(:,:,2)=Y;
        Ic(:,:,3)=I;
    case 6
        Ic(:,:,1)=I;
        Ic(:,:,2)=I;
        Ic(:,:,3)=Y;
    case 7
        Ic(:,:,1)=I;
        Ic(:,:,2)=X;
        Ic(:,:,3)=X;
    case 8
        Ic(:,:,1)=X;
        Ic(:,:,2)=I;
        Ic(:,:,3)=X;
    case 9
        Ic(:,:,1)=X;
        Ic(:,:,2)=X;
        Ic(:,:,3)=I;
        
end
imshow(uint8(Ic))
