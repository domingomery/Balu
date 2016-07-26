% [Az,Sn,Sp1,t] = Bev_roc(z,d,show)
%
% Toolbox: Balu
%    ROC Analysis for feature z with classification c. show = 1 indicates
%    that the ROC curve will be displayed.
%    Az is the area under the ROC curve
%    Sn and Sp1 are the coordinates of Sensitibity and 1-Specificity of 
%    optimal point of ROC curve (nearest point to (0,1).
%    t is the threshold for this optimal point.
%
%  Example:
%     load datagauss          % simulated data (2 classes, 2 features)
%     Bev_roc(X(:,1),d,1)     % ROC plot of first feature of X
%
%     Az = Bev_roc(X(:,2),d); % Az of ROC for second feature of X
%
%  See also Bio_plotfeatures.
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl
%
function [Az,Sn,Sp1,t,x,y,th] = Bev_roc(z,d,show)

if ~exist('show','var')
    show = 0;
end

dmin = min(d);
dmax = max(d);

dn = dmax-dmin+1;

if (dn==2)
    d = d-dmin;


    n = 100;

    x  = ones(n+1,1); %1 - specificity
    y  = ones(n+1,1); %sensibility

    m = length(z);

    ii1 = find(d==1);
    ii0 = find(d==0);

    if (mean(double(z(ii0)))>mean(double(z(ii1))))
        z = -z;
    end


    zmin = min(z);
    zmax = max(z);

    zdelta = (zmax-zmin)/n;
    Az = 0;
    i = 1;

    C0 = zeros(n+3,1);
    C1 = zeros(n+3,1);

    zi = zmax+zdelta:-zdelta:zmin-zdelta;

    for th = zi
        ii = find(z>=th);
        c  = zeros(m,1);
        if (~isempty(ii))
            c(ii) = ones(length(ii),1);
        end
        TP = c'*d;
        TN = (1-c)'*(1-d);
        FP = c'*(1-d);
        FN = (1-c)'*(d);
        x(i,1)  = FP/(FP+TN); %1-specificity
        y(i,1)  = TP/(TP+FN); %sensibility
        if (i>1)
            Az = Az + (x(i,1)-x(i-1,1))*(y(i,1)+y(i-1,1))/2;
            C1(i,1) = sum( (z(ii1)>=th)&(z(ii1)<zant));
            C0(i,1) = sum( (z(ii0)>=th)&(z(ii0)<zant));
        end
        i = i + 1;
        zant = th;
    end

    if (show)
        clf;
        subplot(1,2,1)
        plot(x,y);
        hold on
        plot(x,y,'r.')
        grid on
        text(0.55,0.55,sprintf('Az = %f',Az));
        title('ROC curve')
        xlabel('1-Sp');
        ylabel('Sn');
        axis([-0.1 1.1 -0.1 1.1])
        grid on


        subplot(1,2,2)
        C1n = C1/sum(C1);
        C0n = C0/sum(C0);
        if (Az>0)
            s0 = sprintf('class %d',dmin);
            s1 = sprintf('class %d',dmax);
            plot(zi,C1n,'r.-',zi,C0n,'b'),legend(s1,s0)
            title('PDF')
            xlabel('z');
            ylabel('h(z)')
        end
    end
    
    th = zi;

    % Optimal point of ROC
    r        = sqrt(x.*x + (1-y).*(1-y));
    [rmin,i] = min(r);
    Sn       = y(i);
    Sp1      = x(i);
    t        = zi(i);

else
    disp('I cannot compute ROC for more than two classes :(');
end