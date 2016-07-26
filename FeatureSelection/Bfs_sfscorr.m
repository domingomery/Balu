% [R,selec] = Bfa_sfscorr(f,d,m)
%
% Toolbox: Balu
%    Sequential Forward Selection for features f according to measurment d.
%    The algorithm searchs the linear combination of the features that
%    best correlates with d.
%    m features will be selected.
%    R is the obtained correlation coefficient and selec are the number of the
%    selected features.
%
% D.Mery, PUC-DCC, Aug. 2008
% http://dmery.ing.puc.cl
%

function [R,selec] = Bfa_sfscorr(f,d,m)

selec  = []; %selected features
R    = [];
Rmax = 0;
k=1;


while (k<=m)
    nuevo = 0;
    for i=1:size(f,2)
        if (k==1) || (sum(selec==i)==0)
            if std(f(:,i))>0
                s = [selec i];
                fs = f(:,s);
                [per,Rs] = Bfa_corrsearch(fs,d,1,2,0);
                if (Rs>Rmax)
                    ks = i;
                    Rmax = Rs;
                    nuevo = 1;
                end
            end
        end
    end
    if (nuevo)
        selec = [selec ks];
        R = [R Rmax];
        clf
        bar(R)
        hold on
        for i=1:length(selec)
            text(i-0.4,R(i)*1.05,sprintf('%d',selec(i)));
        end
        pause(0)
        k = k + 1;
    else
        disp('no more improvement, the sequantial search is interrupted.');
        k = 1e10;
    end
end
figure
Bfa_corrsearch(f(:,selec),d,1,2,1);
