% [New_X, New_d] = Bds_TomekLinks(X, d, k)
%
% Toolbox: Balu
%
%    UNDER-SAMPLING method, where the label for the majority class is 0 
%    and the label for the minority class is 1.
%
%    This function implements the Tomek Links Algorithm as an undersampling 
%    method. Only the points belonging to the majority class (labeled 0) in 
%    the links are removed. The links are verified in a neighborhood with 
%    k neighbours. It means that a maximum of k points of the majority 
%    class will be removed by each point of the minority class.
%
%    Tomek link: Ei, Ej belong to different classes, d(Ei, Ej) is the 
%    distance between them. A (Ei, Ej) pair is called a Tomek link if 
%    there is no example El, such that d(Ei, El) < d(Ei, Ej) 
%    or d(Ej , El) < d(Ei, Ej).

%
%    Input:
%       X is a matrix with features (columns).
%       d is the ideal classification for X. 
%       k is the number of neighbours to use in the test
%
%    Output:
%       New_X is the new feature matrix with less elements that X.
%       New_d is the new ideal classification for New_X.
%
%    For more details on the theoretical description of the original 
%    algorithm  please refer to the following paper:
%    I. Tomek (1976). Two modi?cations of CNN. IEEE Transactions on 
%    Systems, Man and Cybernetics, 6(11):769–772, 1976.
%
% C. Mera, UNAL, 2013


function [new_X, new_d] = Bds_TomekLinks(X, d, k)

    % Points of the majority class (labeled 0)
    NEG_DAT = X(d == 0,:);

    % Points of the minority class (labeled 1)
    POS_DAT = X(d == 1,:);
    pos_size = size(POS_DAT, 1);

    % The data is re-order
    X = [POS_DAT;NEG_DAT];
    d = ones(size(d,1),1);
    d(pos_size+1:end,:) = 0;
    
    POS_DAT = POS_DAT';

    % Finding the k nearest neighbours
    I = nearestneighbour(POS_DAT, X', 'NumberOfNeighbours', k+1);
    I = I';

  
    DEL = zeros(size(d,1),1);
    
    for i=1:pos_size
        % Si entre los k vecinos hay puntos de la clase maj, entonces se
        % eliminan aquellos que no los antecede un vecino de la clase min 
        % y que por tanto forman un TOMEK LINK con el punto minoritario i.
        W = I(i,2:end);
        maj = sum(d(W)==0);
        if (maj > 0)
            for l=1:k
                if (d(I(i,l+1)) == 0)
                    DEL(I(i,l+1)) = 1;
                else
                    break;
                end
            end
        end
    end
    
    X((DEL==1),:) = [];
    d((DEL==1),:) = [];

    new_X = X;
    new_d = d;
    
end