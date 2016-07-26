% ds      = Bcl_AdaBoostM1(X,d,Xt,options)  Training & Testing together
% options = Bcl_AdaBoostM1(X,d,options)     Training only
% ds      = Bcl_AdaBoostM1(Xt,options)      Testing only
%
% Toolbox: Balu
%    AdaBoost M1 classifier for two classes
%
%    Design data:
%       X is a matrix with features (columns).
%       d is the ideal classification for X, only of two classes {0,1}.
%       options.iter is the number of iterations.
%       options.weakLearner is a Balu struct of a classifier (LDA by default)
%       options.ClassDist True or false. true indicates that the class
%                    distribution is maintained while doing weighted 
%                    resampling false indicates that the class distribution
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.H is a 1 x iter struct with the weakClassifier training.
%       options.B is a iter x iter vector with Beta values.
%       options.string is a string that describes the performed
%       classification (in this case 'AdaBoostM1').
%
%    Implementation based on the source code of: Barnan Das
%    www.eecs.wsu.edu/~bdas1
%
%    For more details on the theoretical description of the algorithm 
%    please refer to the following paper:
%    Polikar, R. (2006): Ensemble Based Systems in Decision Making. IEEE
%    Circuits ans Systems Magazine, Third Quarter, 21-45.
%
%    Example: Training & Test together:
%       load datagauss                  % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)           % plot feature space
%       op.iter = 10;                   % 10 iterartions
%       ds = Bcl_AdaBoostM1(X,d,Xt,op); % AdaBoost classifier
%       p = Bev_performance(ds,dt)      % performance on test data
%
%    Example: Training only
%       load datagauss                % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)         % plot feature space
%       op.iter = 10;                 % prior probability for each class
%       op.weakLearner.name = 'qda';  % WeakClassifier (QDA BaluClassifier)
%       op.weakLearner.p = [];        % WeakClassifier parameters
%       op = Bcl_AdaBoostM1(X,d,op);  % AdaBoost - training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_AdaBoostM1(Xt,op);   % AdaBoost - testing
%       p = Bev_performance(ds,dt)    % performance on test data
%
% C. Mera, UNAL, 2013
%

function [ds,options] = Bcl_AdaBoostM1(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});

% +CarlosM
% WeakLearner
if ~isfield(options,'weakLearner')
    % LDA by Default
    weakLearner(1).name  = 'lda';    
    weakLearner(1).options.p = [];
else
    weakLearner = options.weakLearner;
end


% +CarlosM 
% Boosting Iterations
if ~isfield(options,'iter')
    T = 10;
else
    T = options.iter;
end

% +CarlosM
if ~isfield(options,'show')
    options.show = 1;
end

% +CarlosM 
%                    true or false. true indicates that the class
%                    distribution is maintained while doing weighted 
%                    resampling and before SMOTE is called at each 
%                    iteration. false indicates that the class distribution
if ~isfield(options,'ClassDist')
    ClassDist = true;
else
    ClassDist = options.ClassDist;
end

options.string = 'AdaBoostM1';

if train
    
    [N,m] = size(X);

    %Positive Class (or Minority Class) == 1
    %Negative Class (or Mayority Class) == 0    
    POS_DATA = X(d == 1,:);
    NEG_DATA = X(d == 0,:);
    
    pos_size = size(POS_DATA, 1);
    neg_size = size(NEG_DATA, 1);

    % Reorganize TRAIN by putting all the positive and negative exampels
    % together, respectively.
    X = [POS_DATA;NEG_DATA];
    d(1:pos_size)   = 1;
    d(pos_size+1:N) = 0;
    
    
    % W stores the weights of the instances in each row for every iteration of
    % boosting. Weights for all the instances are initialized by 1/m for the
    % first iteration.
    W = ones(N,1)/N;
    
    % L stores pseudo loss values, H stores hypothesis, B stores (1/beta) 
    % values that is used as the weight of the % hypothesis while forming the 
    % final hypothesis. % All of the following are of length <=T and stores 
    % values for every iteration of the boosting process.
    L = [];
    H = {};
    B = [];

    % Keeps counts of the number of times the same boosting iteration have been
    % repeated
    count = 0;
    
    % Loop counter
    t = 1;

    % LOG MESSAGE
    if (options.show ~= 0)
        fprintf('AdaBoostM1. WeakClassifier=%s Iteraciones=%d\n', weakLearner.name, T);
        disp (['   ' int2str(t)]);
    end
    
    % Boosting T iterations
    while t <= T

        % LOG MESSAGE
        if (options.show ~= 0)
            disp (['   Boosting iteration #' int2str(t)]);
        end
    
        if ClassDist == true
            % Resampling POS_DATA with weights of positive example
            sum_POS_WT = sum(W(1:pos_size,t));
            POS_WT = W(1:pos_size,t)/sum_POS_WT;
            
            RND_IDX = randsample(1:pos_size,pos_size,true,POS_WT);
            RESAM_POS = POS_DATA(RND_IDX,:);

            % Resampling NEG_DATA with weights of positive example
            sum_NEG_WT = sum(W(pos_size+1:N,t));
            NEG_WT = W(pos_size+1:N,t)/sum_NEG_WT ;

            RND_IDX = randsample(1:neg_size,neg_size,true,NEG_WT);
            RESAM_NEG = NEG_DATA(RND_IDX,:);

            % Resampled TRAIN is stored in RESAMPLED
            RESAMPLED = [RESAM_POS; RESAM_NEG];
            D = zeros(size(RESAM_POS,1)+size(RESAM_NEG,1),1);
            D(1:size(RESAM_POS,1)) = 1;
        else 
            % Indices of resampled train
            RND_IDX = randsample(1:N, N, true, W(:,t));

            % Resampled TRAIN is stored in RESAMPLED
            RESAMPLED = X(RND_IDX,:);
            D = d(RND_IDX,:);

            % Calulating the percentage of boosting the positive class. 'pert'
            % is used as a parameter of SMOTE
            pos_size = sum(RESAMPLED(:,end)==1);
            neg_size = sum(RESAMPLED(:,end)==0);
        end
        
        %Train the Weak Classifier
        model = Bcl_structure(RESAMPLED, D, weakLearner);
        
        %Classify the Training Data with the model obtained to calcule the
        %loos
        pred = Bcl_structure(X, model);
        
        loss = sum(W(d~=pred,t));

        % If count exceeds a pre-defined threshold (5 in the current
        % implementation), the loop is broken and rolled back to the state
        % where loss > 0.5 was not encountered.
        if count > 5
           L = L(1:t-1);
           H = H(1:t-1);
           B = B(1:t-1);
           options.iter = t-1;
           disp ('          Too many iterations have loss > 0.5');
           disp ('          Aborting boosting...');
           break;
        end

        % If the loss is greater than 1/2, it means that an inverted
        % hypothesis would perform better. In such cases, do not take that
        % hypothesis into consideration and repeat the same iteration. 'count'
        % keeps counts of the number of times the same boosting iteration have
        % been repeated
        if loss > 0.5
            count = count + 1;
            continue;
        else
            count = 1;
        end        

        L(t) = loss; % Pseudo-loss at each iteration
        H{t} = model; % Hypothesis function   
        beta = loss/(1-loss); % Setting weight update parameter 'beta'.
        B(t) = log(1/beta); % Weight of the hypothesis

        % At the final iteration there is no need to update the weights any
        % further
        if t==T
            break;
        end

        % Updating weight        
        W(:,t+1) = W(:,t);
        W(d==pred,t+1) = W(d==pred,t)*beta;
        
        % Normalizing the weight for the next iteration
        sum_W = sum(W(:,t+1));
        
        % This verification is necessary in the case that all examples be
        % rigth classified
        if (sum_W > 0)
            W(:,t+1) = W(:,t+1)/sum_W;
        else
            W(:,t+1) = W(:,t);
        end
        
        % Incrementing loop counter
        t = t + 1;
    end    
      
    % Normalizing B
    sum_B = sum(B);
    B = B/sum_B;
    
    %Output Parameters
    options.B = B;
    options.H = H;
    ds = options;
end


if test

    H = options.H;
    B = options.B;
    nt = size(Xt,1); % Total number of instances in the test set
    ht = zeros(nt, T);
    prediction = ones(nt,1);
    for t = 1:T
        ht(:,t) = Bcl_structure(Xt,H{t});
    end

    for i=1:nt
        wt_zero = sum(B(ht(i,:)==0));
        wt_one  = sum(B(ht(i,:)==1));

        if (wt_one > wt_zero)
            prediction(i,:) = 1;
        else
            prediction(i,:) = 0;
        end
    end
    
    ds = prediction;
end


