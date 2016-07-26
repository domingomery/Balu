% [X,Xn] = Bfx_haralick(I,R,options)
% [X,Xn] = Bfx_haralick(I,options)
%
% Toolbox: Balu
%    Haralick texture features.
%
%    X is a 28 elements vector with mean and range of mean and range of
%
%     1 Angular Second Moment
%     2 Contrast
%     3 Correlacion
%     4 Sum of squares
%     5 Inverse Difference Moment
%     6 Sum Average
%     8 Sum Entropy
%     7 Sum Variance
%     9 Entropy
%    10 Difference Variance
%    11 Difference Entropy
%    12,13 Information Measures of Correlation
%    14 Maximal Corrleation Coefficient
%
%    Xn is the list of name features.
%
%    I is the image. R is the binary image that indicates which pixels of I will be
%    computed.
%    options.dharalick is the distance in pixels used to compute the
%    coocurrence matrix.
%    options.show = 1 display results.
%
%    Reference:
%    Haralick (1979): Statistical and Structural Approaches to Texture,
%    Proc. IEEE, 67(5):786-804
%
%   Example 1: only one distance (3 pixels)
%      options.dharalick = 3;                 % 3 pixels distance for coocurrence
%      I = imread('testimg1.jpg');            % input image
%      R = Bim_segbalu(I);                    % segmentation
%      J = I(:,:,2);                          % green channel
%      [X,Xn] = Bfx_haralick(J,R,options);    % Haralick features
%      Bio_printfeatures(X,Xn)
%
%   Example 2: five distances (1,2,...5 pixels)
%      options.dharalick = 1:5;                 % 3 pixels distance for coocurrence
%      I = imread('testimg1.jpg');            % input image
%      R = Bim_segbalu(I);                    % segmentation
%      J = I(:,:,2);                          % green channel
%      [X,Xn] = Bfx_haralick(J,R,options);    % Haralick features
%      Bio_printfeatures(X,Xn)
%
%   See also Bfx_gabor, Bfx_clp, Bfx_fourier, Bfx_dct, Bfx_lbp.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl


function [X,Xn] = Bfx_haralick(I,R,options)

I = double(I);

if nargin==2;
    options = R;
    R = ones(size(I));
end

if isempty(R)
    R = ones(size(I));
end    

dseq = options.dharalick;

if ~isfield(options,'show')
    options.show = 0;
end

if options.show == 1
    disp('--- extracting Haralick texture features...');
end

m = length(dseq);
n = 28*m;

X  = zeros(1,n);
Xn = char(zeros(n,24));
k = 1;
for i=1:m
    d = dseq(i);

    Cd000 = Bcoocurrencematrix(I, R, d, 0)+Bcoocurrencematrix(I,R, -d, 0);Cd000 = Cd000/sum(Cd000(:));
    Cd045 = Bcoocurrencematrix(I, R, d,-d)+Bcoocurrencematrix(I,R, -d, d);Cd045 = Cd045/sum(Cd045(:));
    Cd090 = Bcoocurrencematrix(I, R, 0, d)+Bcoocurrencematrix(I,R,  0,-d);Cd090 = Cd090/sum(Cd090(:));
    Cd135 = Bcoocurrencematrix(I, R, d, d)+Bcoocurrencematrix(I,R, -d,-d);Cd135 = Cd135/sum(Cd135(:));

    TexMat = [Bcoocurrencefeatures(Cd000) Bcoocurrencefeatures(Cd045) Bcoocurrencefeatures(Cd090) Bcoocurrencefeatures(Cd135)];
    X(1,i*28-27:i*28) = [mean(TexMat,2); max(abs(TexMat'))']';

    for q=1:2
        if (q==1)
            sq = 'mean ';
        else
            sq = 'range';
        end
        for s=1:14
            Xn(k,:) = sprintf('Tx%2d,d%2d(%s)         ',s,d,sq);
            k = k + 1;

        end
    end
end
end


% P = Bcoocurrencematrix(I,R,Io,Jo)
%
% Coocurrence matrix of the pixels of image I indicated by binary image R
% following the direction (Io,Jo).
%
% (c) D.Mery, PUC-DCC, Apr. 2008

function P = Bcoocurrencematrix(I,R,Io,Jo)
V = fix(I/32)+1;
[N,M] = size(I);
Z1 = zeros(N+40,M+40);
Z2 = Z1;
R1 = Z1;
R2 = R1;
Z1(15:N+14,15:M+14) = V;
Z2(15+Io:N+14+Io,15+Jo:M+14+Jo) = V;
R1(15:N+14,15:M+14) = R;
R2(15+Io:N+14+Io,15+Jo:M+14+Jo) = R;
ii = find(not(and(R1,R2)));
Z1(ii) = -ones(length(ii),1);
Z2(ii) = -ones(length(ii),1);
T1 = Z1(:);
T2 = Z2(:);
d = find(and((T1>-1),(T2>-1)));
if (not(isempty(d)))
    P = zeros(8,8);
    X = sortrows([T1(d) T2(d)]);
    i1 = find(or(([0; X(:,1)]-[X(:,1); 0]~=0),...
        ([0; X(:,2)]-[X(:,2); 0]~=0)));
    i2 = [i1(2:length(i1)); 0];
    d  = i2-i1;
    for i=1:length(d)-1
        P(X(i1(i),2),X(i1(i),1)) = d(i);
    end
else
    P = -ones(8,8);
end
end


% Tx = Bcoocurrencefeatures(P)
%
% Haralick texture features calculated from coocurrence matrix P.
%
% (c) D.Mery, PUC-DCC, Apr. 2008
function Tx = Bcoocurrencefeatures(P)

Pij = P(:);
Ng = 8;
pxi = sum(P,2);
pyj = sum(P)';

ux = mean(pxi);
uy = mean(pyj);

sx = std(pxi);
sy = std(pyj);


pxy1 = zeros(2*Ng-1,1);
for k=2:2*Ng
    s = 0;
    for i=1:Ng
        for j=1:Ng
            if (i+j == k)
                s = s + P(i,j);
            end
        end
    end
    pxy1(k-1) = s;
end


pxy2 = zeros(Ng,1);
for k=0:Ng-1
    s = 0;
    for i=1:Ng
        for j=1:Ng
            if (abs(i-j) == k)
                s = s + P(i,j);
            end
        end
    end
    pxy2(k+1) = s;
end


Q = zeros(Ng,Ng);
pxi = pxi+1e-20;
pyj = pyj+1e-20;

for i=1:Ng
    for j=1:Ng
        s = 0;
        for k=1:Ng
            s = s + P(i,k)*P(j,k)/pxi(i)/pyj(k);
        end
        Q(i,j) = s;
    end
end
eigQ = eig(Q);


[i,j] = find(P>=0);
dif   = i-j;
dif2 = dif.*dif;
dif21 = dif2 + 1;


% 1 Angular Second Moment
f1 = Pij'*Pij;

% 2 Contrast
f2 = ((0:Ng-1).*(0:Ng-1))*pxy2;

% 3 Correlacion
f3 = (sum(i.*j.*Pij)-ux*uy*Ng^2)/sx/sy;

% 4 Sum of squares
f4 = dif2'*Pij;

% 5 Inverse Difference Moment
f5 = sum(Pij./dif21);

% 6 Sum Average
f6 = (2:2*Ng)*pxy1;

% 8 Sum Entropy
f8 = -pxy1'*log(pxy1+1e-20);

% 7 Sum Variance
if8 = (2:2*Ng)'-f8;
f7  = if8'*pxy1;

% 9 Entropy
f9 = -Pij'*log(Pij+1e-20);

% 10 Difference Variance
f10 = var(pxy2);

% 11 Difference Entropy
f11 = -pxy2'*log(pxy2+1e-20);

% 12,13 Information Measures of Correlation
HXY = f9;
pxipyj = pxi(i).*pyj(j);
HXY1 = -Pij'*log(pxipyj+1e-20);
HXY2 = -pxipyj'*log(pxipyj+1e-20);
HX = -pxi'*log(pxi+1e-20);
HY = -pyj'*log(pyj+1e-20);
f12 = (HXY-HXY1)/max([HX HY]);
f13 = (1-exp(-2*(HXY2-HXY)));

% 14 Maximal Corrleation Coefficient
f14 = (eigQ(2));


Tx = [f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14]';
end
