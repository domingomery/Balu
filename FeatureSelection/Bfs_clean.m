% selec = Bfs_clean(X,show)
%
% Toolbox: Balu
%    Feature selection cleaning.
%   
%    It eliminates constant features and correlated features.
%
%    Input: X is the feature matrix.
%           show = 1 displays results (default show=0)
%    Output: selec is the indices of the selected features
%
% Example:
%    load datareal
%    s = Bfs_clean(f,1);            % index of selected features
%    X = f(:,s);                    % selected features
%    Xn = fn(s,:);                  % list of feature names
%    disp('Original:'); howis(f)    % original set of features
%    disp('Selected:'); howis(X)    % set of selecte features
%
% D.Mery, PUC-DCC, Jul. 2009
% http://dmery.ing.puc.cl


function selec = Bfs_clean(X,show)

f = X;

if not(exist('show','var'))
    show = 0;
end

nf = size(f,2);
p  = 1:nf;
ip = zeros(nf,1);

% eliminating correlated features
warning off
C = abs(corrcoef(f));
warning on
[ii,jj] = find(C>0.99);
if (not(isempty(ii)))
    for i=1:length(ii)
        if (abs(ii(i)-jj(i))>0)
            k = max([ii(i) jj(i)]);
            t = find(p==k);
            n = length(p);
            if (not(isempty(t)))
                if (t==1)
                    p = p(2:n);
                else
                    if (t==n)
                        p = p(1:n-1);
                    else
                        p = p([1:t-1 t+1:n]);
                    end
                end
            end
        end
    end
end

ip(p) = 1;

% eliminating constant features
s = std(f);
ii = find(s<1e-8);
if not(isempty(ii))
    ip(ii) = 0;
end
p = find(ip);
fc  = f(:,p);
nc = size(fc,2);
if show
    fprintf('Bfs_clean: number of features reduced from %d to %d.\n',nf,nc)
end
selec=p;
