% [X,Xn] = Bfx_moments(R,options)
%
% Toolbox: Balu
%
%    Extract moments and central moments.
%
%    options.show    = 1 display mesagges.
%    options.central = 1 for central moments, and 0 for normal moments
%    options.rs      = n x 2 matrix that contains the indices of the
%                      moments to be extracted (see example).
%
%      X is vector that contains the n moments.
%      Xn is the list of the n feature names.
%
%    Example (Centroid of a region)
%      I = imread('testimg1.jpg');     % input image
%      R = Bim_segbalu(I);             % segmentation
%      imshow(R);
%      options.show = 1;
%      options.central = 0;
%      options.rs = [0 0; 1 0; 0 1];
%      [X,Xn] = Bfx_moments(R,options);
%      ic = X(2)/X(1);
%      jc = X(3)/X(1);
%      hold on
%      plot(jc,ic,'rx')
%
%   See also Bfx_basicgeo, Bfx_gupta, Bfx_fitellipse, Bfx_flusser, 
%   Bfx_hugeo.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl

function [X,Xn] = Bfx_moments(R,options)
if ~exist('options','var')
    options.show = 0;
end

if options.show == 1
    disp('--- extracting moments...');
end
[Ireg,Jreg] = find(R==1);           % pixels in the region
A   = length(Ireg);

central = options.central;
rs      = options.rs;
n = size(rs,1);

if central
    i_m = mean(Ireg);
    j_m = mean(Jreg);
    I1 = Ireg - i_m*ones(A,1);
    J1 = Jreg - j_m*ones(A,1);
    strc = 'Central Moment';
else
    I1 = Ireg;
    J1 = Jreg;
    strc = 'Moment';
end

X  = zeros(1,n);
Xn = char(zeros(n,24));

for i=1:n
    r = rs(i,1);
    s = rs(i,2);
    Ir = ones(A,1);
    Js = ones(A,1);
    if r>0
        for jr = 1:r
            Ir = Ir.*I1;
        end
    end
    if s>0
        for js = 1:s
            Js = Js.*J1;
        end
    end
    X(i) = Ir'*Js;
    str = [strc ' ' num2str(r) ',' num2str(s) '                '];
    Xn(i,:) = str(1:24);
end


