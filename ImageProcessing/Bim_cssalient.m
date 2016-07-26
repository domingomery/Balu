% function [J_on, J_off] = Bim_cssalient( I, type, show )
%
% Toolbox: Balu
%   compute center surround saliency image from and gray scale image follow
%   mehtod mentioned in Montabone & Soto (2010). Method precompute an 
%   integral image to do all calculations.
%
%   Inputs:
%   I:  gray scale image. If depth of I is higher that 1 the algoritm
%       transform the image into gray scale version.
%   type:   scales method
%   show:   plot options. 1 shows both saliency images. 0 shows anything.
%
%   Output:
%   J_off: Surround-Center
%   J_on: Center-Surround
%
%   Example:
%   I = imread('testimg7.bmp');
%   [J_on,J_off] = Bim_cssalient(I,1,1);
%
%   Reference:
%   Montabone, S.  and Soto, A. (2010): Human detection using a mobile 
%   platform and novel features derived from a visual saliency mechanism. 
%   Image and Vision Computing. 28(3):391?402.
%
% (c) Christian Pieringer Baeza, cppierin@uc.cl, 2009.

function [J_on, J_off] = Bim_cssalient( I, type, show )

[N, M, D] = size( I );

if D > 1
    I = rgb2gray(I);
end

I = double(I);
SAT = cumsum(cumsum(I), 2);

switch type
    case 1
        % values in Sebastian thesis
        sigma = [12, 24, 28, 48, 56, 112];
        N_sigma = length( sigma );
        
    case 2
        W = max(N,M);
        sigma = [W/8, W/4, W/2];
        N_sigma = length( sigma );
        
%     case 3
%         sigma = [7, 9, 11, 13, 17, 19];
%         N_sigma = length( sigma );
end

J_on = zeros( N, M, N_sigma );
J_off = zeros( N, M, N_sigma );

for scales = 1:N_sigma
    s = sigma(scales)+1;
    SubWindow = zeros(N,M);
    for x = 1:N
        for y = 1:M
            x_min = max( 2, x-s );
            x_max = min( N, x+s );
            y_min = max( 2, y-s );
            y_max = min( M, y+s );
      
            SubWindow(x,y) = SAT( x_min-1, y_min-1 ) ...
            + SAT( x_max, y_max ) ...
            - SAT( x_min-1, y_max ) ...
            - SAT( x_max, y_min-1 );              
        end
    end
    
    Surround = (SubWindow - I) ./ repmat( ( (2*s+1)^2 )-1, N,M);

    J_on(:,:,scales) = max( (I - Surround), 0);
    J_off(:,:,scales) = max( (Surround - I), 0 );

end

% add along through the scales
J_on = sum( J_on, 3 );
J_off = sum( J_off, 3 );

if show
    figure, imshow( I, [] ), title('Input image')
    figure, imshow( J_on, [] ), title('J\_on')
    figure, imshow( J_off, [] ), title('J\_off')
end

end