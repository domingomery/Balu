% [New_X, New_d] = Bds_ROS(X, d)
%
% Toolbox: Balu
%
%    OVER-SAMPLING method, where the label for the majority class is 0 
%    and the label for the minority class is 1.
%
%    This function implments a Random Over Sampling for the minority class 
%   (labeled 1)in a imbalace data set. It takes N samples form minority 
%    class with replacement. N is the number of elements of the majority class.
%
%    Input:
%       X is a matrix with features (columns).
%       d is the ideal classification for X. 
%
%    Output:
%       New_X is the new feature matrix with more elements that X.
%       New_d is the new ideal classification for New_X.
%
%    For more details on the theoretical description of the algorithm 
%    please refer to the following paper:
%    H. He y E. A. Garcia (2009). Learning from Imbalanced Data. IEEE Trans
%    Knowl Data Eng, vol. 21, No.9, pp. 1263–1284, sep. 2009.
%
% C. Mera, UNAL, 2013

function [new_X, new_d] = Bds_ROS(X, d)

    % Points of the majority class (labeled 0)
    NEG_DAT = X(d == 0,:);

    % Points of the minority class (labeled 1)
    POS_DAT = X(d == 1,:);

    neg_size = size(NEG_DAT, 1);
    pos_size = size(POS_DAT, 1);
    
    % Only if the minority class has less elements that the majority
    if (pos_size < neg_size)
        % Select pos_size samples form X with replacement,
        RND_IDX = randsample(1:pos_size, neg_size, true);
        new_X   = [POS_DAT(RND_IDX,:); NEG_DAT];
        new_d   = ones(neg_size*2, 1);
        new_d(neg_size+1:end,:) = 0;
    else
        new_X = X;
        new_d = d;
    end
end