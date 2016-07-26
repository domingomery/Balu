%function Bio_showconfusion(C)
%
% Toolbox: Balu
%    Show confusion matrix C in a color 2D representation.
%    
%
% Example:
%
% % Simulation of a 10x10 confussion matrix
% C = rand(10,10)+2*eye(10);C = C/max(C(:));
%
% Bio_showconfusion(C)
% (c) GRIMA-DCCUC, 2013
% http://grima.ing.puc.cl
%


function Bio_showconfusion(C)

minC = min(C(:));
maxC = max(C(:));


[N,M] = size(C);

if minC<0
    error('There are negative values.');
end

if maxC>100
    error('Maximal value of C should be 100.')
end

if N~=M
    error('Matrix should be square.')
end


if maxC<=1
    C = 100*C;
end

jet100 = imresize(jet,[100 3],'nearest');


clf;
T = 256;
imshow(imresize(round(C),[T T],'nearest'),jet100)
axis off
hold on

sq = 256/N;

for i=1:N
    y = T*i/N-sq/2;
    for j=1:N
        x = T*j/N-0.75*sq;
        text(x,y,num2str(round(C(i,j))))
    end
end
colorbar
        