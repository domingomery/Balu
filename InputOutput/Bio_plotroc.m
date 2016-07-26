% Bio_plotroc(FPR,TPR,col)
%
% Toolbox: Balu
%   Plot ROC curve and fit to an exponetial curve
%
% Example:
%
% th = 3;x = 0:0.05:1; y = 1-exp(-3*x)+randn(1,21)*0.05;
% Bio_plotroc(x,y)
%
% D.Mery, PUC-DCC, Apr. 2013
% http://dmery.ing.puc.cl
%

function [AUC,TPRs,FPRs,TPR05] = Bio_plotroc(x,y,col_line,col_point)

if ~exist('col_line','var')
    col_line = 'b';
end

if ~exist('col_point','var')
    col_point = 'r.';
end

% clf
plot(x,y,col_point)

ths = fminsearch(@thest,1,[],x,y);
xs = 0:0.005:1;
a = 1/(1-exp(-ths));
ys = a*(1-exp(-ths*xs));

AUC = a*(1-exp(-ths)/ths-1/ths);

hold on
plot(xs,ys,col_line)

d = xs.*xs + (1-ys).*(1-ys);

[~,ii] = min(d);

TPRs = ys(ii(1));
FPRs = xs(ii(1));

ii = find(xs==0.05);
TPR05 = ys(ii(1));

plot(FPRs,TPRs,[col_point(1) '*'])

xlabel('FPR')
ylabel('TPR')
end



function err = thest(th,x,y)
a  = 1/(1-exp(-th));
ys = a*(1-exp(-th*x));
err = norm(y-ys);
end




