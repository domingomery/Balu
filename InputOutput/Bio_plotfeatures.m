% Bio_plotfeatures(X,d,Xn)
%
% Toolbox: Balu
%    Plot features X acording classification d. If the feature names are
%    given in Xn then they will labeled in each axis.
%    
%    For only one feature, histograms are ploted.
%    For two (or three) features, plots in 2D (or 3D) are given.
%    For m>3 features, m x m 2D plots are given (feature i vs. feature j)
%
% Example 1: 1D & 2D
%    load datagauss                    % simulated data (2 classes, 2 features)
%    figure(1)
%    Bio_plotfeatures(X(:,1),d,'x1')   % histogram of feature 1
%    figure(2)
%    Bio_plotfeatures(X,d)             % plot feature space in 2D (2 features)
%
% Example 2: 3D
%    load datareal                     % real data 
%    X = f(:,[221 175 235]);           % only three features are choosen
%    Bio_plotfeatures(X,d)             % plot feature space in 3D (3 features)
%
% Example 3: 5D (using feature selection)
%    load datareal                     % real data
%    op.m = 5;                         % 5 features will be selected
%    op.s = 0.75;                      % only 75% of sample will be used
%    op.show = 0;                      % display results
%    op.b.name = 'fisher';             % definition SFS with Fisher
%    s = Bfs_balu(f,d,op);             % feature selection      
%    Bio_plotfeatures(f(:,s),d)        % plot feature space for 5 features
%
% See also Bev_roc.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl
%

function Bio_plotfeatures(X,d,Xn)

if ~exist('Xn','var')
    Xn = [];
end

m = size(X,2);

if m>9
    error('Bio_plotfeatures for %d features makes %d plots.',m,m*m)
end

scflag = 0;
if size(d,2)==2
    sc = d(:,2);
    d  = d(:,1);
    scflag = 1;
end


dmin = min(d);
dmax = max(d);




col = 'gbrcmykbgrcmykbgrcmykbgrcmykgbrcmykbgrcmykbgrcmykbgrcmykgbrcmykbgrcmykbgrcmykbgrcmyk';
mar = 'ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^';
% clf
warning off
r = sprintf('%s',39);
warning on
s = 'legend(';
for k = dmin:dmax
    s = [s r sprintf('class %d',k) r ];
%    s = [s r sprintf('class %d',k-1) r ];
    if k<dmax
        s = [s ','];
    else
        s = [s ');'];
    end
end
if (m<4)
    switch m
        case 1
            for k = dmin:dmax
                [h,x] = hist(X(d==k),100);
                dx = x(2)-x(1);
                A = sum(h)*dx;
                x = [x(1)-dx x x(end)+dx];
                h = [0 h 0];
                plot(x,h/A,col(k+1));
                hold on
            end
            if ~isempty(Xn)
                xlabel(Xn);
            else
                xlabel('feature value');
            end
            ylabel('PDF');
        case 2
            if isempty(Xn)
                Xn = ['feature value 1';'feature value 2'];
            end
            for k = dmin:dmax
                ii = find(d==k);
                plot(X(ii,1),X(ii,2),[col(k+1) mar(k+1)]);
                hold on
                if scflag
                    for ic=1:length(ii)
                        text(X(ii(ic),1),X(ii(ic),2),['  ' num2str(sc(ii(ic)))]);
                    end
                end
            end
            title('feature space');
            xlabel(Xn(1,:));
            ylabel(Xn(2,:));
        case 3
            for k = dmin:dmax
                ii = find(d==k);
                plot3(X(ii,1),X(ii,2),X(ii,3),[col(k+1) mar(k+1)]);
                hold on
                if scflag
                    for ic=1:length(ii)
                        text(X(ii(ic),1),X(ii(ic),2),X(ii(ic),3),['  ' num2str(sc(ii(ic)))]);
                    end
                end
            end
            if ~isempty(Xn)
                xlabel(Xn(1,:));
                ylabel(Xn(2,:));
                zlabel(Xn(3,:));
            else
                xlabel('feature value 1');
                ylabel('feature value 2');
                zlabel('feature value 3');
            end
    end
    eval(s)
else
    l = 1;
    for j=1:m
        for i=1:m
            zi = X(:,i);
            zj = X(:,j);
            subplot(m,m,l); l = l+1;

            for k = dmin:dmax
                ii = find(d==k);
                plot(zi(ii),zj(ii),[col(k+1) mar(k+1)]);
                hold on
                if scflag
                    for ic=1:length(ii)
                        text(zi(ii(ic)),zj(ii(ic)),['  ' num2str(sc(ii(ic)))]);
                    end
                end
                
            end
            if ~isempty(Xn)
                xl = Xn(i,:);
                yl = Xn(j,:);
            else
                xl = sprintf('z_%d',i);
                yl = sprintf('z_%d',j);
            end
            if i==1
                ylabel(yl)
            end
            if j==m
                xlabel(xl)
            end
        end
    end


end