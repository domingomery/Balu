function Bio_plotsquare(varargin)
% (c) Domingo Mery and German Mondragon
n = nargin;
c = varargin{n};
if compare(class(c),'char')==0
    n = n-1;
else
    c = 'b';
end

switch n
    case 1
        ix = varargin{1};
    case 2
        ix1 = varargin{1};
        ix2 = varargin{2};
        ix = [ix1(1) ix2(1) ix1(2) ix2(2)];
    case 4
        ix = [varargin{1} varargin{2} varargin{3} varargin{4}];
    otherwise
        error('Bio_plotsquare does not work with this input');
end
hold on
if (size(ix,1)==1)
    i1 = ix(1);
    i2 = ix(2);
    j1 = ix(3);
    j2 = ix(4);
    x = [j1 j2 j2 j1 j1];
    y = [i1 i1 i2 i2 i1];
    plot(x,y,c)
    % plot(x,y,'g','LineWidth',2)
else
    for ind=1:size(ix,1)
        i1 = ix(ind, 1);
        i2 = ix(ind, 2);
        j1 = ix(ind, 3);
        j2 = ix(ind, 4);
        x = [j1 j2 j2 j1 j1];
        y = [i1 i1 i2 i2 i1];
        plot(x,y,c)
        % plot(x,y,'g','LineWidth',2)
    end
end
