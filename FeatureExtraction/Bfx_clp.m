% [X,Xn] = Bfx_clp(I,R,options)
% [X,Xn] = Bfx_clp(I,options)
%
% Toolbox Balu: Crossing Line Profile.
%
%    X is the features vector, Xn is the list of feature names(see Example
%    to see how it works).
%
%    Reference:
%    Mery, D.: Crossing line profile: a new approach to detecting defects
%    in aluminium castings. Proceedings of the Scandinavian Conference on
%    Image Analysis 2003 (SCIA 2003), Lecture Notes in Computer Science
%    LNCS 2749: 725-732, 2003.
%
%   Example:
%      options.show    = 1;                    % display results
%      options.ng      = 32;                   % windows resize
%      I = imread('testimg4.jpg');             % input image
%      J = I(395:425,415:442,1);               % region of interest (red)
%      R = J>135;                              % segmentation
%      figure;imshow(J,[])
%      figure;imshow(R)
%      [X,Xn] = Bfx_clp(J,R,options);          % CLP features
%      Bio_printfeatures(X,Xn)
%
%   See also Xcontrast, Xharalick, Bfx_clp, Xfourier, Xdct, Xlbp.
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [X,Xn] = Bfx_clp(I,R,options)

[N,M] = size(I);
I = double(I);

if nargin==2;
    options = R;
    R = ones(N,M);
end

if ~isfield(options,'show')
    options.show = 0;
end

if options.show == 1
    disp('--- extracting Crossing line profile features...');
end

show = options.show;



[ii,jj] = find(R==1);
h  = max(ii)-min(ii)+1;                    % height
w  = max(jj)-min(jj)+1;                    % width
x  = round(Xcentroid(R));                  % mass center
i1 = max([1 x(1)-h]); j1 = max([1 x(2)-w]);
i2 = min([N x(1)+h]); j2 = min([M x(2)+w]);

%if ~isempty(R);
%    I(R==0) = 0;
%end

ng = options.ng;

Bn = imresize(I(i1:i2,j1:j2),[ng ng]);

mg = fix(ng/2+1);

% Crossing line profiles
P0 = Bn(mg,:)';     %   0.0 ...  90.0  4
P1 = Bn(:,mg);      %  90.0 ...   0.0  0
P2 = zeros(ng,1);   %  45.0 ... 135.0  6
P3 = zeros(ng,1);   % 135.0 ...  45.0  2
P4 = zeros(ng,1);   %  22.5 ... 112.5  5
P5 = zeros(ng,1);   %  67.5 ... 157.5  7
P6 = zeros(ng,1);   % 112.5 ...  22.5  1
P7 = zeros(ng,1);   % 157.5 ...  67.5  3


Q0 = Bn; Q1 = Bn; Q2 = Bn; Q3 = Bn;
Q4 = Bn; Q5 = Bn; Q6 = Bn; Q7 = Bn;

Q0(mg,:) = 255*ones(1,ng);
Q1(:,mg) = 255*ones(ng,1);

m4 = mg/(ng-1);
b4 = mg/2-m4;
b7 = 3*mg/2+mg/(ng-1);


for i=1:ng
    P2(i,1) = Bn(i,i);
    Q2(i,i) = 255;
    
    P3(i,1) = Bn(i,ng-i+1);
    Q3(i,ng-i+1) = 255;
    
    j4      = fix(m4*i + b4 + 0.5);
    P4(i,1) = Bn(j4,i);
    Q4(j4,i) = 255;
    
    P5(i,1) = Bn(i,j4);
    Q5(i,j4)  = 255;
    
    j7      = fix(-m4*i + b7 + 0.5);
    P6(i,1) = Bn(i,j7);
    Q6(i,j7)  = 255;
    
    P7(i,1) = Bn(j7,i);
    Q7(j7,i) = 255;
end


PP = [P0 P1 P2 P3 P4 P5 P6 P7];

d = abs(PP(1,:)-PP(ng,:));

[~,J] = sort(d);

Po = PP(:,J(1));

Po = Po/Po(1);

m  = (Po(ng)-Po(1))/(ng - 1);

mb = Po(1)-m;

Q = Po-(1:ng)'*m-ones(ng,1)*mb;

Qm = mean(Q);

Qd  = max(Q)-min(Q);

Qd1 = log(Qd+1);
Qd2 = 2*Qd/(Po(1)+Po(ng));
Qs = std(Q);

Qf = fft(Q);
Qf = abs(Qf(2:8,1));

if (show)
    
    figure(10)
    clf
    subplot(2,4,5);plot(PP(:,1));axis([1 ng 0 255]);title('k=4');
    subplot(2,4,1);plot(PP(:,2));axis([1 ng 0 255]);title('k=0');
    subplot(2,4,3);plot(PP(:,3));axis([1 ng 0 255]);title('k=2');
    subplot(2,4,7);plot(PP(:,4));axis([1 ng 0 255]);title('k=6');
    subplot(2,4,4);plot(PP(:,5));axis([1 ng 0 255]);title('k=3');
    subplot(2,4,2);plot(PP(:,6));axis([1 ng 0 255]);title('k=1');
    subplot(2,4,8);plot(PP(:,7));axis([1 ng 0 255]);title('k=7');
    subplot(2,4,6);plot(PP(:,8));axis([1 ng 0 255]);title('k=5');
    figure(11)
    imshow([Q1 Q5 Q2 Q4;Q0 Q7 Q3 Q6],gray(256));
    pause(0);
end

X = [Qm Qs Qd Qd1 Qd2 Qf'];

Xn = [ 'CLP-Qm                  '
    'CLP-Qs                  '
    'CLP-Qd                  '
    'CLP-Qd1                 '
    'CLP-Qd2                 '
    'CLP-Qf1                 '
    'CLP-Qf2                 '
    'CLP-Qf3                 '
    'CLP-Qf4                 '
    'CLP-Qf5                 '
    'CLP-Qf6                 '
    'CLP-Qf7                 '];

