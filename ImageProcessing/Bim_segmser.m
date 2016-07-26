% [frames,seeds,M] = Bim_segmser(I,param)
%
% Toolbox: Balu
%  Segmentation using MSER algorithm.
%
%     I: input image
%     param(1) = minimal area of the ellipse in pixels
%     param(2) = maximal area of the ellipse in pixels
%     param(3) = MinDiversity eg 0.7
%     param(4) = MaxVariation eg 0.2
%     param(5) = Delta eg 10
%     param(6) = 0: DarkOnBright, 1: BrightOnDark; 2: Both
%     param(7) = 1: show results
%
%     frames is the output according vl_mser, i.e., each column corresponds
%     to the center (x0,y0), the radius and the orientation of the
%     segmented ellipse.
%
%    M is a matrix with the same size as I whose value are equal to the
%    number of overlapping extremal regions. It is like a segmented image.
% 
% Reference: J. Matas, O. Chum, M. Urban, and T. Pajdla, "Robust wide
%       baseline stereo from maximally stable extremal regions," in
%       Proc. BMVC, 2002.
% 
%  Example:
%     I  = imread('X1.png');
%     [fr,sd,J] = Bim_segmser(I,[5 2500 0.7 0.2 10 0 1]);
%     figure
%     imshow(J,[])
%
%  See also vl_mser, vl_erfill.
%
% (c) GRIMA-DCC-PUC, 2011
% http://dmery.ing.puc.cl

function [frames,seeds,M] = Bim_segmser(I,param)

Amin = param(1); % minimal area of the ellipse
Amax = param(2); % maximal area of the ellipse
Dmin = param(3); % MinDiversity eg 0.7
Vmax = param(4); % MaxVariation eg 0.2
Pdel = param(5); % Delta eg 10
B_D  = param(6); % 0: DarkOnBright, 1: BrightOnDark; 2: Both
show = param(7);


if show
    clf
    imshow(I,[]);
end

switch B_D
    case 0 % Dark objects only
        [r,f1] = vl_mser(uint8(I),'MinDiversity',Dmin,'MaxVariation',Vmax,'Delta',Pdel,'BrightOnDark',0) ;
    case 1 % Bright objects only
        [r,f1] = vl_mser(uint8(I),'MinDiversity',Dmin,'MaxVariation',Vmax,'Delta',Pdel,'DarkOnBright',0) ;
    case 2 % both
        [r,f1] = vl_mser(uint8(I),'MinDiversity',Dmin,'MaxVariation',Vmax,'Delta',Pdel,'BrightOnDark',1,'DarkOnBright',1) ;
end

f2 = vl_ertr(f1) ;
fr = [];
rr = [];
for i=1:size(f2,2)
    xo = f2(1,i);
    yo = f2(2,i);
    s1 = f2(3,i);
    s2 = f2(4,i);
    s3 = f2(5,i);
    A = s1==0;
    B = s2==0;
    C = s3==0;
    if not(and(B,or(A,C)))

        % Ellipse equation: X'*inv(S)*X = 1, con X=[x-xo;y-yo] S = [s1 s2;s2 s3]
        % (s3*x^2)/(s1*s3 - s2^2) + (s3*xo^2)/(s1*s3 - s2^2) + (s1*y^2)/(s1*s3 - s2^2) + (s1*yo^2)/(s1*s3 - s2^2) - (2*s3*x*xo)/(s1*s3 - s2^2) -(2*s2*x*y)/(s1*s3 - s2^2) + (2*s2*x*yo)/(s1*s3 - s2^2) + (2*s2*xo*y)/(s1*s3 - s2^2) - (2*s2*xo*yo)/(s1*s3 - s2^2) - (2*s1*y*yo)/(s1*s3 - s2^2) - 1 = 0
        % transformation a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
        % a = s3/(s1*s3-s2^2);
        % b = -(2*s2)/(s1*s3-s2^2);
        % c = s1/(s1*s3-s2^2);
        % d = -(2*s3*xo)/(s1*s3-s2^2)+(2*s2*yo)/(s1*s3-s2^2);
        % e = (2*s2*xo)/(s1*s3-s2^2)-(2*s1*yo)/(s1*s3-s2^2);
        % f = (s3*xo^2)/(s1*s3 - s2^2)+(s1*yo^2)/(s1*s3-s2^2)-(2*s2*xo*yo)/(s1*s3-s2^2)-1;

        % transformation (a*x^2 + b*x*y + c*y^2 + d*x + e*y + f)/(s1*s3-s2^2) = 0
        a = s3;
        b = -(2*s2);
        c = s1;
        d = -(2*s3*xo)+(2*s2*yo);
        e = (2*s2*xo)-(2*s1*yo);
        f = (s3*xo^2)+(s1*yo^2)-(2*s2*xo*yo)-(s1*s3-s2^2);

        aa = [a b c d e f];

        a = aa/aa(6);

        % get ellipse orientation
        alpha = atan2(a(2),a(1)-a(3))/2;

        % get scaled major/minor axes
        ct = cos(alpha);
        st = sin(alpha);
        ap = a(1)*ct*ct + a(2)*ct*st + a(3)*st*st;
        cp = a(1)*st*st - a(2)*ct*st + a(3)*ct*ct;

        % get translations
        T = [[a(1) a(2)/2]' [a(2)/2 a(3)]'];
        mc = -inv(2*T)*[a(4) a(5)]';

        % get scale factor
        val = mc'*T*mc;
        scale = abs(1 / (val- a(6)));

        % get major/minor axis radii
        ae  = 1/sqrt(scale*abs(ap));
        be  = 1/sqrt(scale*abs(cp));
        % ecc = ae/be; % eccentricity (not used)
        ar  = pi*ae*be;

 
        if (ar>=Amin)&&(ar<=Amax)
            fr = [fr [xo yo sqrt(ar) alpha-pi/2]'];
            rr = [rr; r(i)];
            if show
                vl_plotframe(f2(:,i)) ;
                enterpause(0)
            end
        end
    end
end
frames = fr;
seeds = rr;
M = zeros(size(I)) ;
Iu = uint8(I);
for x=rr'
 s = vl_erfill(Iu,x) ;
 M(s) = M(s) + 1;
end
