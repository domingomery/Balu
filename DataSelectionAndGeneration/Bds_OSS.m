% [New_X, New_d] = Bds_OSS(X, d)
%
% Toolbox: Balu
%
%    UNDER-SAMPLING method, where the label for the majority class is 0 
%    and the label for the minority class is 1.
%
%    This function implements the ONE SIDE SELECTION algorithm to 
%    undersampling the majority class (labeled 0). First the CNN Rule is
%    aplied and then all examples of the negative class that participate on
%    a Tomek Link are removed.
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
%    M. Kubat and S. Matwin (1997). Addressing the curse of imbalanced 
%    training sets: One-sided selection. In Proc. 14th International 
%    Conference on Machine Learning, pages 179–186, 1997.
%
% C. Mera, UNAL, 2013


function [new_X, new_d] = Bds_OSS(X, d)

    % First, Use CNN Rule
    [A, B] = Bds_CNNRule(X,d);
    
    % The number max number of neg examples to erase by each pos data is
    % calculated. It isn't part of the original algorithm but neccesary for
    % the implementation of Tomek Links
    pos_size = sum(B == 1);
    neg_size = sum(B == 0);
    
    k = round(neg_size/pos_size);
    
    if (k < 1)
        k =1;
    end
    
    % Then, erase all negative sample on a Tomek Link
    [new_X, new_d] = Bds_TomekLinks(A, B, k);
    
end