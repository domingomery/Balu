% function [param]=fitSigmoid2Scores(posFeats, negFeats, classifier, display)
% original de A.Soto
% nPosExamples=size(posFeats,1);
% nNegExamples=size(negFeats,1);
%
% [predict_label, accuracy, posScores]=svmpredict31EchoOff(ones(nPosExamples,1), posFeats, classifier);
% [predict_label, accuracy, negScores]=svmpredict31EchoOff(ones(nNegExamples,1), negFeats, classifier);
%
% Modificado por D.Mery
%
function [ysc,param] = Bft_sigmoid(sc,d,show)

if ~exist('show','var')
    show = 0;
end


if length(d)>2 % searching for fit

    dneg         = min(d);
    dpos         = max(d);

    posScores    = sc(d==dpos);
    negScores    = sc(d==dneg);

    nPosExamples = length(posScores);
    nNegExamples = length(negScores);

    x            = [negScores; posScores];

    y            = [zeros(nNegExamples,1); ones(nPosExamples,1)];

    param        = nlinfit(x,y,@fsig,[-1 0.05]);
else
    x            = sc;
    param        = d;
end

ysc              = fsig(param,x);


if show
    figure;
    plot(x,y,'bo');
    hold;
    plot(posScores,fsig(param,posScores),'g*');
    plot(negScores,fsig(param,negScores),'r*');
    [i,j] = sort(ysc);
    plot(x(j),ysc(j),'k')
end;


end


function y = fsig(param,x)

y = 1./(1 + exp(param(1)*x+param(2)));

end

%param=fit(x,y,'1./ (1 + exp(a*x+b))','start',[0 20]);