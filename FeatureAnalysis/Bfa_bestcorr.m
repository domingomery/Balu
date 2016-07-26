% [selec,cx,a] = Bfa_bestcorr(X,y)
%
% Toolbox: Balu
%    Search of the variables of X that best correlate with y.
%
%    selec is the number of the selected variables
%    cx   is the correlation coefficient
%    a is the parameters vector of the model
%    selec and cx are sorted (eg, selec(i) is the number of the
%    column of X that corresponds to the i-th best variable
%    that correlates with y, and cx(i) is the corresponding
%    correlation coefficient.
%
%    if correlation with features^2 is requered it can be used
%    the following command: [selec,cx] = Bfa_bestcorr([X X.*X],y);
%
%    Example: 
%       N = 100;
%       x1 = rand(N,1);              % correlated with x1
%       x2 = rand(N,1);              % not correlated with x1
%       x3 = x1 + 0.01*rand(N,1);    % correlated with x1
%       x4 = x1 + 0.001*rand(N,1);   % correlated with x1v
%       X = [x1 x2 x3 x4];
%       y  = x1 + 0.0001*rand(N,1);  % correlated with x1
%       [selec,cx,a] = Bfa_bestcorr(X,y)
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl
%
function [selec,cx,a] = Bfa_bestcorr(X,y)
close all

m = size(X,2);
cc = zeros(m,1);
for i=1:m
    x = X(:,i);
    if std(x)>0
        cc(i) = abs(corr(x,y));
    end
end

[cx,selec] = sort(cc,'descend');

x = X(:,selec(1));
plot(x,y,'.');
hold on
xlabel('best correlated feature (x1)');
ylabel('measured value (y)');

n = length(x);
XX = [x ones(n,1)];

a = (XX'*XX)\XX'*y;   % a = inv(XX'*XX)*XX'*y;

s1 = sprintf('y=%f*x1+%f',a(1),a(2));


ax = axis;
XX = [ax(1:2)' ones(2,1)];
YY = XX*a;
hold on
plot(XX(:,1),YY,'r')

title([sprintf('abs(R)=%5.4f, ',cx(1)) s1])

