% [X,Xn] = Bfx_all(I,R,options)
% [X,Xn] = Bfx_all(I,options)
%
% Toolbox: Balu
%    All pixels.
%
%    X is the features vector, Xn is the list feature names (see Example to
%    see how it works).
%
%
%   Example:
%      I = imread('testimg1.jpg');             % input image
%      [X,Xn] = Bfx_all(I);                    % all pixels
%
%
% (c) D.Mery, PUC-DCC, 2012
% http://dmery.ing.puc.cl

function [X,Xn] = Bfx_all(I,R,options)


X = I(:)';

Xn = zeros(length(X),24);
