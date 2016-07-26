%  L = Bim_regionseed(I,seeds)
%
% Toolbox: Balu
%  Growing region algorithm from seed pixels.
%  I    : grayscale image
%  seeds: seed pixels
%
%  Example:
%     I = imread('X1.png');
%     [fr,seeds] = Bim_segmser(I,[5 2500 0.7 0.2 10 0 1]);
%     L = Bim_regionseed(I,seeds);
%     figure;imshow(L,[])
%
% D.Mery, PUC-DCC, Dec. 2011
% http://dmery.ing.puc.cl
%

function L = Bim_regionseed(I,seeds,N)
I = double(I);
[NI,MI] = size(I);
L = zeros(NI,MI);
Dt = L;
if ~exist('N','var')
    N = 8;
end
if N~=8
    X = [0 1 0; 1 1 1; 0 1 0];
else
    X = ones(3,3);
end
for r=1:length(seeds)
    D = zeros(NI,MI);
    [i,j] = ind2sub([NI MI],seeds(r));
    m = I(i,j);
    %J = and(I<=(1.25*m),not(Dt));
    J = and(I<=m,not(Dt));
    % clf;
    % imshow(J,[]);
    % hold on
    % plot(j,i,'rx');
    % enterpause
    D(i,j) = 1;
    ok = 0;
    nd0 = 1;
    while not(ok)
        nd = nd0;
        D = and(J,imdilate(D,X));
        nd0 = sum2(D);
        ok = nd0-nd==0;
    end
    Dt = or(D,Dt);
    L = L + double(D)*r;
end
