% [New_X, New_d] = Bds_NCL(X, d)
%
% Toolbox: Balu
%
%    UNDER-SAMPLING method, where the label for the majority class is 0 
%    and the label for the minority class is 1.
%
%    This function implements the  Neighborhood  Cleaning  Rule to 
%    undersampling the majority class (labeled 0). The majority class  
%    examples that invade minority class space and lead to misclassification
%    of minority class examples in that region are removed. 
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
%    J. Laurikkala (2001). Improving identi?cation of difficult small 
%    classes by balancing class distribution. In Conference on arti?cial 
%    intelligence in medicine in Europe No 8, Portugal, Vol 2101, 63–66, 2001.
%
% C. Mera, UNAL, 2013


function [new_X, new_d] = Bds_NCL(X, d)

    % Finding the 3 nearest neighbours of all samples
    I = nearestneighbour(X', X', 'NumberOfNeighbours', 4);
    I = I';
    
    INDX = zeros(size(d));
    
    for i=1 : size(d,1)
        W = I(i,2:end);
        p = sum(d(W)==1);
        
        % If i belongs to the majority class and the classification given 
        % by its 3-NN contradicts the original class of i, then i is removed.
        if (d(i)==0 && p >= 2)
            INDX(i) = 1;
        end

        % If i belongs to the minority class and its 3-NN misclassify i, 
        % then the NNs that belong to the majority class are removed
        if (d(i)==1 && p < 2)
            INDX(W(d(W)==0)) = 1;
        end
    end
    
    X(INDX==1,:) = [];
    d(INDX==1,:) = [];
    
    new_X = X;
    new_d = d;
end