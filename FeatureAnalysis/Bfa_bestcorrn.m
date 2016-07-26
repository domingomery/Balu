% [selec,cx,B] = Bfa_bestcorrn(X,y)
%
% Toolbox: Balu
%    Search of the variables of X that best correlate with y.
%
%    Three models are obtained:
%          y  = a1*z1 + a0 > see figure 1
%          y  = a1*z1 + a2*z2 + a0 > see figure 2
%          y  = a1*z1 + a2*z2 + a3*z3 + a0 > see figure 3
%
%    selec is the number of the selected variables, i.e., 
%       z1 = X(:,selec(1))
%       z2 = X(:,selec(2))
%       z3 = X(:,selec(3))
%    for model i, aj is stored in B(i,j+1) 
%    cx is the correlation coefficient for each model
%    selec and cx are sorted (eg, selec(i) is the number of the
%    column of X that corresponds to the i-th best variable
%    that correlates with y, and cx(i) is the corresponding
%    correlation coefficient.
%
%    if correlation with features^2 is requered it can be used
%    the following command: [selec,c] = Bfa_bestcorrn([X X.*X],y);
%
%       N = 100;
%       x1 = rand(N,1);                    % correlated with x1
%       x2 = rand(N,1);                    % not correlated with x1
%       x3 = x1 + 0.01*rand(N,1);          % correlated with x1
%       x4 = x1 + 0.001*rand(N,1);         % correlated with x1v
%       X = [x1 x2 x3 x4];
%       y  = x1 + 3*x2 + 0.001*rand(N,1);  % correlated with x1
%       [selec,cx,B] = Bfa_bestcorrn(X,y)
%
% D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl
%

function [selec,cx,B] = Bfa_bestcorrn(X,y)
close all

B = zeros(4,3);

p = y;

[N,M] = size(X);

% 1) Best correlated features: Model y = a1*z1 + a0
[selec,c,a] = Bfa_bestcorr(X,p);
k1 = selec(1);
cx(1) = c(1);
B(1,1:2) = a(1:2);

% Best correlated features: Model y = a1*z1 + a2*z2 + a0
c2 = zeros(M,1);
for k2= 1:M
    if (k1~=k2)
        Xb = [X(:,[k1 k2]) ones(N,1)];
        b = regress(p,Xb);
        ps = Xb*b;
        c2(k2)=abs(corr(ps,p));
    end
end
[i,j] = max(c2);
cx(2) = i;
k2 = j;


figure(2)
Xb = [X(:,[k1 k2]) ones(N,1)];
b = regress(p,Xb);
ps = Xb*b;
x = ps;
y = p;
plot(x,y,'.');
hold on
hold on
xlabel('best correlated features (z1+z2+1)');
ylabel('measured value (y)');
s1 = sprintf('y=%f*z1+%f*z2+%f',b(1),b(2),b(3));
B(2,1:3) = b(1:3);

n = length(x);
XX = [x ones(n,1)];
a = (XX'*XX)\XX'*y; % a = inv(XX'*XX)*XX'*y;

ax = axis;
XX = [ax(1:2)' ones(2,1)];
YY = XX*a;
hold on
plot(XX(:,1),YY,'r')


title(sprintf('abs(R)=%5.4f, %s',cx(2), s1))


% Best correlated features: Model y = a1*z1 + a2*z2 + a3*z3 + a0

c3 = zeros(M,1);

for k3= 1:M

    if and((k3~=k2),(k3~=k1))
        Xb = [X(:,[k1 k2 k3]) ones(N,1)];
        b = regress(p,Xb);
        ps = Xb*b;
        c3(k3)=abs(corr(ps,p));
    end
end
[i,j] = max(c3);
cx(3) = i;

k3 = j;
selec = [k1 k2 k3];

Xb = [X(:,[k1 k2 k3]) ones(N,1)];
b = regress(p,Xb);
ps = Xb*b;


figure(3)
clf
x = ps;
y = p;
plot(x,y,'.');
hold on
xlabel('best correlated features (z1+z2+z3+1)');
ylabel('measured value (y)');
s1 = sprintf('y=%f*z1+%f*z2+%f*z3+%f',b(1),b(2),b(3),b(4));
B(3,1:4) = b(1:4);

n = length(x);
XX = [x ones(n,1)];
a = (XX'*XX)\XX'*y; % a = inv(XX'*XX)*XX'*y;

ax = axis;
XX = [ax(1:2)' ones(2,1)];
YY = XX*a;
hold on
plot(XX(:,1),YY,'r')

title(sprintf('abs(R)=%5.4f, %s',i,s1))

