% [p,ci] = Bev_crossval(X,d,options)
%
% Toolbox: Balu
%
%    Cross-validation evaluation of a classifier.
%
%    v-fold Cross Validation in v groups of samples X and classification d
%    according to given method. If v is equal to the number of samples,
%    i.e., v = size(X,1), this method works as the original
%    cross-validation, where training will be in X without sample i and
%    testing in sample i. ci is the confidence interval.
%
%    X is a matrix with features (columns)
%    d is the ideal classification for X
%
%    options.b is a Balu classifier or several classifiers (see example)
%    options.v is the number of groups (folders) of the cross-validations
%    options.c is the probability of the confidence intervale.
%    options.strat = 1 means the portions are stratified.
%    options.show = 1 displays results.
%
%    Example for one classifier without stratified data:
%       load datagauss                                                        % simulated data (2 classes, 2 features)
%       Bplotfeatures(X,d)                                                    % plot feature space
%       b.name = 'knn'; b.options.k = 5;                                      % knn with 5 neighbors
%       op.strat = 0; op.b = b; op.v = 10; op.c = 0.90; op.show = 0;     	  % 10 groups cross-validation for 90%
%       [p,ci] = Bev_crossval(X,d,op)                                         % cross valitadion
%
%    Example for more classifiers with stratified data:
%       load datagauss                                                        % simulated data (2 classes, 2 features)
%       b(1).name = 'knn';   b(1).options.k = 5;                              % KNN with 5 neighbors
%       b(2).name = 'knn';   b(2).options.k = 7;                              % KNN with 7 neighbors
%       b(3).name = 'knn';   b(3).options.k = 9;                              % KNN with 9 neighbors
%       b(4).name = 'lda';   b(4).options.p = [];                             % LDA
%       b(5).name = 'qda';   b(5).options.p = [];                             % QDA
%       b(6).name = 'nnglm'; b(6).options.method = 3; b(6).options.iter = 10; % Nueral network
%       b(7).name = 'libsvm';b(7).options.kernel = '-t 2';                    % rbf-SVM
%       b(8).name = 'dmin';  b(8).options = [];                               % Euclidean distance
%       b(9).name = 'maha';  b(9).options = [];                               % Mahalanobis distance
%       op.strat=1; op.b = b; op.v = 10; op.show = 1; op.c = 0.95;        	  % 10 groups cross-validation
%       [p,ci] = Bev_crossval(X,d,op);                                        % cross valitadion
%
%    The mean performance of classifier k is given in p(k), and the
%    confidence intervals for c*100% are in ci(k,:).
%
% (c) Grima, PUC-DCC, 2011: D.Mery, E.Cortazar, B. Banerjee
% http://grima.ing.puc.cl

function [p,ci] = Bev_crossval(X,d,options)

v        = options.v;
b        = options.b;
show     = options.show;
c        = options.c;

if isfield(options,'strat')
    strat = options.strat;
else
    strat = 0;
end

if (v==1)
    disp('Warning: cross validation with only one group means data training = data test.');
end

if not(exist('show','var'))
    show=1;
end


n = length(b);
N = size(X,1);

dmin = min(d);
dmax = max(d);
nn   = dmin:dmax;

p   = zeros(n,1);
ci  = zeros(n,2);
for k=1:n

    if (v==1)
        XX = X;
        XXt = X;
        dd = d;
        ddt = d;
    elseif (strat==1)
		temporal = cell(v,2);
		for cl=dmin:dmax
			selec = (d==cl);
			XTemp = X(selec,:);
			dTemp = d(selec,:);
			
			numElem = sum(selec);
            
            % Checking if every class has at least v samples, otherwise raise Warning (Sandipan)
            if numElem < v
                warning('Class %d has only %d samples, less than fold value %d!!!\n',dTemp(1,1),numElem,v)
            end
            
			[i,j] = sort(rand(numElem, 1));
			r = floor(numElem/v);
			
			XTemp = XTemp(j,:);
			dTemp = dTemp(j);
			
			for iTemp=1:v
				if iTemp == v
					rango = ((iTemp-1)*r + 1):numElem;
				else
					rango = ((iTemp-1)*r + 1):((iTemp)*r);
				end
				temporal{iTemp,1} = [temporal{iTemp,1};XTemp(rango,:)];
				temporal{iTemp,2} = [temporal{iTemp,2};dTemp(rango)];
			end
		end
		Xr = [];
		dr = [];
		R = zeros(v,2);
		ant = 0;
		for iTemp=1:v
			Xr = [Xr;temporal{iTemp,1}];
			dr = [dr;temporal{iTemp,2}];
			next = size(dr,1);
			R(iTemp,:) = [ant+1 next];
			ant = next;
		end
	else
        rn = rand(N,1);
        [i,j] = sort(rn);

        Xr = X(j,:);
        dr = d(j);

        r = fix(N/v);
        R = zeros(v,2);
        ini = 1;
        for i=1:v-1
            R(i,:) = [ini ini+r-1];
            ini = ini + r;
        end
        R(v,:) = [ini N];
    end
	
    pp = zeros(v,1);
    for i=1:v
        if (v>1)
            XXt = Xr(R(i,1):R(i,2),:);
            ddt = dr(R(i,1):R(i,2),:);
            XX = [];
            dd = [];
            for j=1:v
                if (j~=i)
                    XX = [XX;Xr(R(j,1):R(j,2),:)];
                    dd = [dd;dr(R(j,1):R(j,2),:)];
                end
            end
        end
        [dds,ops] = Bcl_structure(XX,dd,XXt,b(k));
        pp(i) = Bev_performance(ddt,dds,nn);
    end

    p(k) = mean(pp);
    s = ops.options.string;
    % Confidence Interval

    if (v>1)
        pm     = mean(pp);
        mu    = pm;
        sigma = sqrt(pm*(1-pm)/N);
        t = (1-c)/2;
        if v>20
            z = norminv(1-t);
        else
            z = tinv(1-t,v-1);
        end
        p1 = max(0,mu - z*sigma);
        p2 = min(1,mu + z*sigma);
        ci(k,:) = [p1 p2];

        if show
            fprintf('%3d) %s  %5.2f%% in (%5.2f, %5.2f%%) with CI=%2.0f%% \n',k,s,p(k)*100,p1*100,p2*100,c*100);
        end
    else
        ci(k,:) = [0 0];
    end
end