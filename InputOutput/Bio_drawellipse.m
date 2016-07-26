% Bio_drawellipse(v,ecol)
%
% Toolbox: Balu
%    Draws an ellipse with a(1)x^2 + a(2)xy + a(3)y^2 + a(4)x + a(5)y + a(6) = 0
%    ecol is the color of the ellipse
%
%    Extracted from http://homepages.inf.ed.ac.uk/rbf/CVonline/
%    CVonline: The Evolving, Distributed, Non-Proprietary, On-Line Compendium
%    of Computer Vision
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl

function Bio_drawellipse(v,ecol)

%  convert to standard form: ((x-cx)/r1)^2 + ((y-cy)/r2)^2 = 1
%   rotated by theta
%   v = solveellipse(a);

% draw ellipse with N points
if not(exist('ecol','var'))
    ecol = 'b';
end

n = size(v,1);
hold on

for r=1:n
    ae    = v(r,3);
    be    = v(r,4);
    theta = v(r,5);
    mcx   = v(r,1);
    mcy   = v(r,2);

    if n>1
            text(mcx,mcy,num2str(r))
    end

    N = 100;
    dx = 2*pi/N;
    R = [ [ cos(theta) sin(theta)]', [-sin(theta) cos(theta)]'];
    X = zeros(N+1);
    Y = X;
    for i = 1:N
        ang = i*dx;
        x = ae*cos(ang);
        y = be*sin(ang);
        d1 = R*[x y]';
        X(i) = d1(1) + mcx;
        Y(i) = d1(2) + mcy;
    end
    X(N+1) = X(1);
    Y(N+1) = Y(1);
    plot(X,Y,ecol)

end