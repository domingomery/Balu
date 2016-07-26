% [M,emean,estd,ecorr] = Bim_labparam(RGBmes,LABmes,model,show)
%
% Toolbox: Balu
%    Estimate the parameters of the conversion RGB->L*a*b after
%    Leon,K.; Mery,D.; Pedreschi,F.;Leon,J.(2006): Color measurement in
%    L*a*b* units from RGB digital images. Food Research International,
%    39(10):1084-1091.
%  
%    RGBmes matrix nx3 with RGB components of n samples
%    LABmes matrix nx3 with L*a*b* components of n samples
%    model 1 is linear 2 is quadratic
%    show = 1 plot the results
%
%    M parameters of the model (see paper)
%    emean, estd error in mean value, standard deviation
%    ecorr is correlation factor
%
% D.Mery, PUC-DCC, May. 2008
% http://dmery.ing.puc.cl
%

function [M,emean,estd,ecorr] = Bim_labparam(RGBmes,LABmes,model,show)
if not(exist('show','var'))
    show = 0;
end

R = RGBmes(:,1);
G = RGBmes(:,2);
B = RGBmes(:,3);
n = size(R,1);
switch model
    case 1 % linear model
        X = [R G B ones(n,1)];
        s = 'linear';
    case 2
        X = [R G B R.*G R.*B G.*B R.*R G.*G B.*B ones(n,1)];
        s = 'quad';
end

M      = [];
emean  = zeros(3,1);
estd   = zeros(3,1);
ecorr  = zeros(3,1);
comp  = ['L*';'a*';'b*'];
for i=1:3
    Y        = LABmes(:,i);
    th       = (X'*X)\X'*Y; % inv(X'*X)*X'*Y;
    M = [M;th'];
    Ys       = X*th;
    R        = corrcoef(Y,Ys);
    emean(i) = mean(abs(Ys-Y));
    estd(i)  = std(abs(Ys-Y));
    ecorr(i) = R(1,2);
    if show
        figure
        plot(Y,Ys,'.')
        xlabel([comp(i,:) '-real']);
        ylabel([comp(i,:) '-modeled ' s]);
        hold on
        ax = axis;
        x = [Y ones(n,1)];
        th = (x'*x)\x'*Ys; % inv(x'*x)*x'*Ys;
        xx = [ax(1) 1;ax(2) 1];
        yy = xx*th;
        hold on
        plot(xx(:,1),yy,'r:')
        T = sprintf('R=%f mean=%f std=%f',ecorr(i),emean(i),estd(i));
        title(comp(i,:))
        text(ax(1)+0.1*(ax(2)-ax(1)),ax(4)-0.2*(ax(4)-ax(3)),T)
    end
end
save LAB RGBmes LABmes model M emean estd ecorr
