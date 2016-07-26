% Bex_arrows
%
% Toolbox: Balu
%    Example: Separation between three types of arrows.
%
%    A very good separability is achieved using moments 2 and 4.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl

clt
% Definitions


% Training
I     = imread('arrows_1.bmp');
J     = imdilate(not(I),ones(3,3));
L     = bwlabel(J,4);

% Results obteined manually with d = Bio_labelregion(I,L,3);
d     = [ones(17,1); 2*ones(12,1); 3*ones(14,1)];

% Feature extraction
b(1).name = 'basicgeo';    b(1).options.show=1;  % basic geometric fetaures
b(2).name = 'hugeo';       b(2).options.show=1;  % Hu moments
b(3).name = 'flusser';     b(3).options.show=1;  % Flusser moments
op.b        = b;

[X,Xn] = Bfx_geo(L,op);                          % geometric features
[X2,Xn2] = Bfs_noposition(X,Xn);                 % delete position features

% Selected features (manually in this example!)

k  = 17; % First Hu Moment
x  = X2(:,k);
xn = Xn2(k,:);

Bio_plotfeatures(x,d,xn)
enterpause

% The first moment of Hu shows a high separability, thus only the Hu
% moments must be extracted for the test.


% Classifier design
m1 = mean(x(d==1));
m2 = mean(x(d==2));
m3 = mean(x(d==3));

th1 = (m1+m2)/2;
th2 = (m2+m3)/2;

% Testing
I     = imread('arrows_2.bmp');
J     = imdilate(not(I),ones(3,3));
[L,n] = bwlabel(J,4);
op.b  = b(2);                                    % only Hu moments

[X,Xn] = Bfx_geo(L,op);                          % geometric features

% Selected features
x = X(:,1);

A1 = zeros(size(I));
A2 = zeros(size(I));
A3 = zeros(size(I));

for i=1:n
    R = L==i;
    if x(i) > th1
        A1 = or(A1,R);
    else
        if x(i)<th2
            A3 = or(A3,R);
        else
            A2 = or(A2,R);
        end
    end
end
figure(2)
imshow(A1); title('Arrows 1');

figure(3)
imshow(A2); title('Arrows 2');

figure(4)
imshow(A3); title('Arrows 3');



