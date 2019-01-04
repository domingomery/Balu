% Bio_manyplot(x,y,labels,xlab,showper)
%
% Toolbox: Balu
%
% Plot of many (x,y) graphs
%
%      Input: 
%       - x could be a column vector with n elements or a matrix with
%         nxm elements
%       - y is a matrix with nxm elements
%      Output:
%      - a plot of m curves, each curve has n points with different colors
%
%
% Example 1: same x for different y
%   clf
%   x  = (0:0.01:1)';
%   y = [sin(x*2*pi) cos(x*pi) sin(x*pi)];
%   Bio_manyplot(x,y)
%   grid on
%   legend({'curve-1','curve-2','curve-3'})
%
% Example 2: different x for different y
%   clf
%   x  = [sort(rand(100,1)) sort(rand(100,1)) sort(rand(100,1))];
%   y = [sin(x(:,1)*2*pi) cos(x(:,2)*pi) sin(x(:,3)*pi)];
%   Bio_manyplot(x,y)
%   grid on
%   legend({'curve-1','curve-2','curve-3'})
%
% (c) Domingo Mery, 2019
% http://dmery.ing.puc.cl

function Bio_manyplot(x,y)

hold on
set(gca,'FontSize',12)
sc = 'bgrcmykbgrcmykbkbgrcmy';
sq = 'sdvo^>sdvo^>vo^>sd';
sl = '-:-:-:-:-:-:-:-:-:-:-:-:-:';
p = size(y,2);
m = size(x,2);
if m==1
    x = repmat(x,[1 p]);
end
for i=1:p
    plot(x(:,i),y(:,i),[sq(i) sl(i)],'LineWidth',2,...
        'MarkerEdgeColor',sc(i),...
        'MarkerFaceColor',sc(i+1),...
        'MarkerSize',5)
end

