% ds      = Bcl_EasyEnsemble(X,d,Xt,options)  Training & Testing together
% options = Bcl_EasyEnsemble(X,d,options)     Training only
% ds      = Bcl_EasyEnsemble(Xt,options)      Testing only
%
% Toolbox: Balu
%
%    EasyEnsemble classifier for imabalance data, where the label for 
%    the majority class is 0 and the label for the minority class is 1.
%
%    It exploit the majority class examples, avoiding 
%    important majority class examples to be ignored by common 
%    under-sampling while maintaining the fast training speed of 
%    under-sampling.  It construct an ensemble of ensembles, where each 
%    individual classi?er is also an ensemble of AdaBoostM1.
%
%
%    Design data:
%       X is a matrix with features (columns).
%       d is the ideal classification for X only of two classes {0,1} 
%         (majority and minority classes, respectively)
%       options.iter is the number of iterations of EasyEnsemble
%       options.roundsAdaBoost is the number of rounds for each AdaBoost
%       options.weakLearner is a Balu structur of a classifier
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.H is a 1 x iter struct with the weakClassifier training.
%       options.B is a iter x iter vector with Beta values.
%       options.string is a string that describes the performed
%       classification (in this case 'BalanceCascade').
%
%    For more details on the theoretical description of the algorithm 
%    please refer to the following paper:
%    Liu, X.-Y. et al (2009) Exploratory Undersampling for Class-Imbalance Learning
%    IEEE Transactions on Systems, Man and Cybernetics - Part B. 39(2): 539-550
%
% C. Mera, UNAL, 2013


function [ds,options] = Bcl_EasyEnsemble(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});

% WeakLearner of AdaBoost
if ~isfield(options,'weakLearner')
    % LDA by Default
    weakLearner(1).name  = 'lda';    
    weakLearner(1).options.p = [];
else
    weakLearner = options.weakLearner;
end

% EasyEnsemble Iterations
if ~isfield(options,'iter')
    T = 10;
else
    T = options.iter;
end

% Boosting Iterations
if ~isfield(options,'roundsAdaBoost')
    options.roundsAdaBoost = 10;
end


% Show a Log?
if ~isfield(options,'show')
    options.show = false;
end

options.string = 'EasyEnsemble';

if train

    %Positive Class (or Minority Class) == 1
    %Negative Class (or Mayority Class) == 0
    POS_DATA = X(d == 1,:);
    NEG_DATA = X(d == 0,:);
    
    pos_size = size(POS_DATA, 1);
    neg_size = size(NEG_DATA, 1);

    % Target vector of resampled data
    D = zeros(2*pos_size,1);
    D(1:pos_size) = 1;

    %WeakCLassifier of AdaBoost
    op.show        = options.show;
    op.iter        = options.roundsAdaBoost;
    op.ClassDist   = true;
    op.weakLearner = weakLearner;
    
    %Vector of AdaBoost ensembles
    H = cell(T);
    
    %Loop counter of EasyEnsemble Iterations
    t = 1;

    while t <= T

        % LOG MESSAGE
        if (options.show ~= 0)
            disp (['EasyEsemble #' int2str(t)]);
        end
    
        % Resampling with replacement of NEG_DATA
        RND_IDX = randsample(1:neg_size, pos_size, true);
        RESAM_NEG = NEG_DATA(RND_IDX,:);

        % Resampled TRAIN is stored in RESAMPLED
        RESAMPLED = [POS_DATA; RESAM_NEG];

        % Train the ensemble with the WeakClassifier
        model = Bcl_AdaBoostM1(RESAMPLED, D, op);
        
        % Store the ensemble Hi
        H{t} = model;
        
        % Incrementing loop counter
        t = t + 1;
    end    
      
    %Output Parameters
    options.H = H;
    ds = options;
end


if test

    H  = options.H;  % Ensemble of AdaBoost Ensembles
    nt = size(Xt,1); % Total number of instances in the test set

    prediction = ones(nt,1);
    HT = [];
    B  = [];
    for t=1:T
            h = H{t}.H;
            b = H{t}.B;
            l = H{t}.iter;
            
            ht = zeros(nt, l);
            
            for j=1:l
                ht(:,j) = Bcl_structure(Xt,h{j});
            end
            
            HT = [HT ht];
            B  = [B b];
    end    
    
    for i=1:nt
        wt_zero = sum(B(HT(i,:)==0));
        wt_one  = sum(B(HT(i,:)==1));
        
        if (wt_one <= wt_zero)
            prediction(i,:) = 0;
        end
    end
    
    ds = prediction;
    
end


