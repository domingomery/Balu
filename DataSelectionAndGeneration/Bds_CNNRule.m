% [New_X, New_d] = Bds_CNNRule(X, d)
%
% Toolbox: Balu
%
%    UNDER-SAMPLING method, where the label for the majority class is 0 
%    and the label for the minority class is 1.
%
%    This function implements the Condensed Nearest Neighbor Rule (CNN) to
%    undersampling the majority class (labeled 0). CNN Rule  aims to ?nd 
%    examples far from decision boundaries. The idea is to ?nd a consistent 
%    data subset Z ? X, where all the examples in X can be classi?ed 
%    correctly by using 1-NN in Z.
%
%    Input:
%       X is a matrix with features (columns).
%       d is the ideal classification for X. 
%
%    Output:
%       New_X is the new feature matrix with less elements that X.
%       New_d is the new ideal classification for New_X.
%
%    For more details on the theoretical description of the original 
%    algorithm  please refer to the following paper:
%    P. E. Hart (1968). The Condensed Nearest Neighbor Rule. 
%    IEEE Transactions on Information Theory IT-14 (1968), 515-516
%
% C. Mera, UNAL, 2013


function [new_X, new_d] = Bds_CNNRule(X, d)

    % Points of the majority class (labeled 0)
    NEG_DAT = X(d == 0,:);
    neg_size = size(NEG_DAT, 1);
    
    % Points of the minority class (labeled 1)
    POS_DAT = X(d == 1,:);
    pos_size = size(POS_DAT, 1);

    % Select one data sample at random from the majority class
    x = randi(neg_size);
    
    % Create Z as all positive samples and one randomly selected negative example
    Z  = [POS_DAT;NEG_DAT(x,:)];
    dz = ones(pos_size+1,1);
    dz(end) = 0;
    
    % Classify X with the 1-NN rule using the examples in Z
    op.k = 1;
    op = Bcl_knn(Z,dz,op);
    
    ds = Bcl_knn(X,op);
    
    INDX = d~=ds;
    Z = [Z; X(INDX,:)];

    new_X = Z;
    new_d = [dz; d(INDX)];
    
end