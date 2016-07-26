% [New_X, New_d] = Bds_smote(X, d, smote_rate, k)
%
% Toolbox: Balu
%
%    OVER-SAMPLING method, where the label for the majority class is 0 
%    and the label for the minority class is 1.
%
%    This function implements SMOTE (Synthetic Minority Over-sampling TEchnique)
%    an oversampling technique for imbalanced data sets, with two classes. 
%    The algorithm generates a synthetic sample  as a point along the line 
%    between each minority sample and one random point selected of their 
%    K-NN in the same class.
%
%    Input:
%       X is a matrix with features (columns).
%       d is the ideal classification for X. 
%       smote_rate is the rate of the synthetic samples to be created for 
%       the minority class.  A value of 1.0 means that 100% of the data in 
%       the minority class will be generated.
%       k is the number of neighbours to use
%
%    Output:
%       New_X is the new feature matrix augmented with the synthetic samples
%       New_d is the new ideal classification for New_X.
%
%    For more details on the theoretical description of the algorithm 
%    please refer to the following paper:
%    Chawla, N. V. et al (2002). SMOTE: synthetic minority over-sampling
%    technique in J. Arti?cial Intell. Res., vol. 16, n.o 1, pp. 321?357
%
%    This implementation use the nearestneighbour provided by Richard
%    Brown,available on: 
%    http://www.mathworks.com/matlabcentral/fileexchange/38830-smote-synthe
%    tic-minority-over-sampling-technique
%
% C. Mera, UNAL, 2013

function [New_X, New_d] = Bds_smote(X, d, smote_rate, k)

    if (nargin < 3)
        smote_rate =   1.0;
        k = 5;
    elseif (nargin < 4)
        k = 5;
    end
    
    if (smote_rate > 0)

        % POS_DAT = candidate points of the minority class (labeled 1)
        POS_DAT = X(d == 1,:);
        POS_DAT = POS_DAT';
        
        pos_size = size(POS_DAT, 2);

        % Calculate the number N of synthetic data examples that need 
        % to be generated for the minority class
        N  = round(pos_size * smote_rate);

        % If the rate es less than 1, only a random sample of the minority
        % class is used
        if (smote_rate < 1)
            RND_IDX = randsample(pos_size, N);
            POS_DAT = POS_DAT(:,RND_IDX);
            pos_size = size(POS_DAT, 2);
        end

        % Finding the k nearest neighbours the positive points
        I = nearestneighbour(POS_DAT, POS_DAT, 'NumberOfNeighbours', k+1);
        I = I';

        [r, c] = size(I);

        Xn = ones(N,1);

        XX = zeros(size(POS_DAT,1), N);
        j = 1;

        % If the rate is a multiple of 1
        while (N >= pos_size)
            for i = 1:r
                index = I(i,randi(k)+1);
                alpha = rand;
                new_P = (1-alpha).*POS_DAT(:,i) + alpha.*(POS_DAT(:,index));
                XX(:,j) = new_P;
                N = N-1;
                j = j+1;
            end
        end

        % The remaning N synthetic samples to be created
        while (N > 0)
            i = randi(r);
            index = I(i, randi(k)+1);
            alpha = rand;
            new_P = (1-alpha).*POS_DAT(:,i) + alpha.*(POS_DAT(:,index));
            XX(:,j) = new_P;
            N = N-1;    
            j = j+1;
        end

        New_X = [X; XX'];
        New_d = [d; Xn];
    else
        New_X = X;
        New_d = d;
    end
end