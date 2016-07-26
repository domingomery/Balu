% Bex_humoments
%
% Toolbox: Balu
%    Example: Separation between 1, 2 and 3 using Hu moments.
%
%    A very good separability is achieved using moments 2 and 4.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl

clt
% Definitions
k           = [2 4];                      % selected features

Itrain      = imread('123.bmp');
Rtrain      = imdilate(not(Itrain),ones(3,3));
[Ltrain,n]  = bwlabel(Rtrain,4);
Ltrain      = bwlabel(Rtrain,8);
imshow(Ltrain,[])
d = [1 1 1 2 1 3 2 2 2 1 2 3 3 3 3]';    % ideal classification
b.name      = 'hugeo';      
b.options.show = 1;                      % basic geometric features
op.b        = b;
[Xtrain,Xn] = Bfx_geo(Ltrain,op);        % geometric features

% Selected features
x   = Xtrain(:,k);
xn  = Xn(k,:);
figure(2)
Bio_plotfeatures(x,d,xn)
