% function H = Bim_inthist(I,b)
%
%
% Toolbox: Balu
%    Integral histogram.
%
%    Input data:
%       I grayvalue image (only for positive values)
%       b number of bins
%
%    Output:
%       H Integral histogram (size = NxMxb, where [N,M] = size(I))
%
%    Example: 
%       I = imread('rice.png');
%       H = Bim_inthist(I,256);
%       i1=70;j1=50;i2=130;j2=190;
%       K = I(i1:i2,j1:j2);
%       t = hist(K(:),1:256);
%       h = Bim_inthistread(H,i1,j1,i2,j2);
%       compare(h,t)
%
%    See also Bim_inthistread.
%
% (c) D.Mery, PUC-DCC, 2012
% http://dmery.ing.puc.cl
% 

function H = Bim_inthist(I,b)
[N,M] = size(I);
J = round(I);
s = log(N*M)/log(2);
if  s > 16
    H = single(zeros(N,M,b));
else
    if s > 8
       H = uint16(zeros(N,M,b));
    else
       H = uint8(zeros(N,M,b));
    end
end
B = J(1,1);
H(1,1,B) = 1;
i = 1;
for j=2:M
    B = J(1,j);
    H(i,j,:) = H(i,j-1,:);
    H(i,j,B) = H(i,j,B) + 1;
end
for i=2:N
    B = J(i,1);
    H(i,1,:) = H(i-1,1,:);
    H(i,1,B) = H(i,1,B) + 1;
    for j=2:M
        B = J(i,j);
        H(i,j,:) = H(i,j-1,:)+H(i-1,j,:)-H(i-1,j-1,:);
        H(i,j,B) = H(i,j,B) + 1;
    end
end



