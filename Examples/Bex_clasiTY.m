% Bex_clasiTY
%
% Toolbox: Balu
%    Example: Separation between T and Y
%
%    This example uses the eccentricity to separate T and Y.
%
%    The classification is performed using a threshold computed as the
%    average of the class centroids.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl

clt

% Definitions
ths = 190;                               % Segmentation threshold
k   = 16;                                % selected feature

% Training
figure(1)
Itrain = imread('TY_1.jpg');
Rtrain = Itrain<ths;                     % segmentation
Ltrain = bwlabel(Rtrain,8);
imshow(Ltrain,[])
d = [1 2 2 1 2 1 2 2 1 1 2]';            % ideal classification
                                         % '1' is Y and '2' is T
b.name         = 'basicgeo';      
b.options.show = 1;                      % basic geometric features
op.b      = b;
[Xtrain,Xn] = Bfx_geo(Ltrain,op);        % geometric features
figure(2)

% Selected feature
x  = Xtrain(:,k);
xn  = Xn(k,:);
Bio_plotfeatures(x,d,xn)
m1  = mean(x(d==1));
m2  = mean(x(d==2));
thc = (m1+m2)/2;                         % classification threshold


% Testing
Itest = imread('TY_2.jpg');
Rtest = Itest<ths;
[Ltest,mtest] = bwlabel(Rtest,8);
[Xtest,Xn] = Bfx_geo(Ltest,op);     % geometric features

ii = Xtest(:,k)>thc;

Y = zeros(size(Ltest));
T = Y;
for i=1:mtest
    if ii(i)==1
        Y = or(Y,Ltest==i);
    else
        T = or(T,Ltest==i);
    end
end

figure(3)
imshow(Y)
figure(4)
imshow(T)