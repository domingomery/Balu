% [New_X, New_d] = Bds_Adasyn(X, d, smote_rate, k)
%
% Toolbox: Balu
%
%    OVER-SAMPLING method, where the label for the majority class is 0 
%    and the label for the minority class is 1.
%
%    This function implements ADASYN (Adaptive Synthetic Sampling Approach) 
%    for Imbalanced Learning. Its main idea is to generate minority class 
%    examples adaptively according to their distributions, where more 
%    synthetic data is created for "difficult" ones than "easy" ones.
%
%    Input:
%       X is a matrix with features (columns).
%       d is the ideal classification for X  {0, 1} where 1 is the
%       minority class.
%       smote_rate is the rate of the synthetic samples to be created for 
%       the minority class.  A value of 1.0 means that 100% of the data in 
%       the minority class will be generated
%       k is the number of neighbours to use
%
%    Output:
%       New_X is the new feature matrix augmented with the synthetic samples
%       New_d is the new ideal classification for New_X.
%
%    For more details on the theoretical description of the algorithm 
%    please refer to the following paper:
%    HE, H. et al (2008). ADASYN: Adaptive synthetic sampling approach for
%    imbalanced learning in IEEE International Joint Conference on NN, 2008
%
%    This implementation use the nearestneighbour provided by Richard Brown
%    and available on: 
%    http://www.mathworks.com/matlabcentral/fileexchange/38830-smote-synthe
%    tic-minority-over-sampling-technique

% C. Mera, UNAL, 2013

function [new_X, new_d] = Bds_Adasyn(X, d, smote_rate, k)

    if (nargin < 3)
        smote_rate =   1.0;
        k = 5;
    elseif (nargin < 4)
        k = 5;
    end
   
    if (smote_rate > 0)
        % Points of the majority class (labeled 0)
        NEG_DAT = X(d == 0,:);
        NEG_DAT = NEG_DAT';

        % Points of the minority class (labeled 1)
        POS_DAT = X(d == 1,:);
        POS_DAT = POS_DAT';

        neg_size = size(NEG_DAT, 2);
        pos_size = size(POS_DAT, 2);

        if (pos_size < neg_size)

            % Calculate the number N of synthetic data examples that need 
            % to be generated for the minority class
            N = round(pos_size * smote_rate);

            % Finding the k nearest neighbours of all the positive samples
            I = nearestneighbour(POS_DAT, X', 'NumberOfNeighbours', k+1);
            I = I';

            % Calculate the ratio r_i 
            R = sum(d(I(:,2:end))==0, 2);
            R = R/k;

            % Normalize r_i
            R = R/sum(R);

            % Calculate the number of synthetic data examples that need to be 
            % generated for each minority example
            G = round(R*N);

            % Coorrect N by the effect of the round function
            N = sum(G(:,1));
            
            XX = zeros(size(POS_DAT,1), N);

            aux = 0;
            for i=1:pos_size
                if (G(i) >= 1)
                    I2 = nearestneighbour(POS_DAT(:,i), POS_DAT, 'NumberOfNeighbours', k+1);
                    for j=1:G(i)
                        index = I2(1,randi(k)+1);
                        aux = aux + 1;
                        alpha = rand;
                        new_P = (1-alpha).*POS_DAT(:,i) + alpha.*(POS_DAT(:,index));
                        XX(:,aux) = new_P;
                    end
                end
            end

            %Delete No-Generated samples
            if (aux < N)
                XX(:,aux+1:end)=[];
                N = aux;
            end

            new_X = [X; XX'];
            new_d = [d; ones(N,1)];

        else
            new_X = X;
            new_d = d;
        end
    else
        new_X = X;
        new_d = d;
    end
end