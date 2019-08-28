% J = Bmv_projective2D(I,H,SJ,show)
%
% Toolbox: Balu
%
%    2D proyective transformation.
%
%    J = projective2D(I,H,SJ,show) returns a new image J that is computed
%    from the 2D projective transformation H of I.
%
%    SJ is [NJ MJ] the size of the transformed image J. The
%    default of SJ is [NJ,MJ] = size(I).
%
%    show = 1 displays images I and J.
%
%    The coordinates of I are (xp,yp) (from x',y'), and the 
%    coordinates of J are (x,y). The relation between both
%    coordinate systems is: mp = H*m, where mp = [xp yp 1]'
%    and m = [x y 1]'.
%
%    R. Hartley and A. Zisserman. Multiple View Geometry in Computer
%    Vision. Cambridge University Press, 2000.
%
%    Example:
%       I = rgb2gray(imread('testimg4.jpg'));  % original image     
%       s = 7.5; t = pi/4;                     % scale and orientation
%       H = [s*cos(t) -s*sin(t)  -500
%            s*sin(t)  s*cos(t) -1750
%               0         0        1];         % similarity matrix
%       J = Bmv_projective2D(I,H,[600 250],1); % projective transformation
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [J,R] = Bmv_projective2D(I,H,SJ,show)

I = double(I);

[NI,MI] = size(I);

if ~exist('SJ','var')
    [NJ,MJ] = size(I);
else
    NJ = SJ(1);
    MJ = SJ(2);
end

if ~exist('show','var')
    show = 0;
end


if (show)
    figure(1)
    imshow(I,[])
    title('image original')
    axis on
end


J = mean(I(:))*ones(NJ,MJ);
R = zeros(NJ,MJ);

X = (1:NJ)'*ones(1,MJ);
Y = ones(NJ,1)*(1:MJ);

x = X(:);
y = Y(:);
m = [x'; y'; ones(1,NJ*MJ)];
mp = H*m;
mp = mp./(ones(3,1)*mp(3,:));
mp = fix(mp + [0.5 0.5 0]'*ones(1,NJ*MJ));
mm = [mp(1:2,:); x'; y'];
mm = mm(:,mm(1,:)>0);
mm = mm(:,mm(2,:)>0);
mm = mm(:,mm(1,:)<=NI);
mm = mm(:,mm(2,:)<=MI);
xp = mm(1,:);
yp = mm(2,:);
x  = mm(3,:);
y  = mm(4,:);

i  = xp + (yp-1)*NI;
j  = x  + (y-1)*NJ;

J(j) = I(i);
R(j) = 1;

if (show)
    figure(2)
    imshow(J,[])
    axis on
    title('transformed image')
end

