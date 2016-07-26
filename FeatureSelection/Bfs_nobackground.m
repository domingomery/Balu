% [f_new,fn_new] = Bfs_nobackground(f,fn)
%
% Toolbox: Balu
%     This procedure deletes the features related to the position.
%     It deletes the following features related to the contrast:
%     - contrast-K1             
%     - contrast-K2             
%     - contrast-K3             
%     - contrast-Ks             
%     - contrast-K              
%
% Example:
%        options.show    = 1;                    % display results
%        options.neighbor = 2;                   % neigborhood is imdilate
%        options.param    = 5;                   % with 5x5 mask
%        I = imread('testimg4.jpg');             % input image
%        J = I(395:425,415:442,1);               % region of interest (red)
%        R = J<=130;                             % segmentation
%        figure;imshow(J,[])
%        figure;imshow(R)
%        [X1,Xn1] = Bfx_contrast(J,R,options);   % contrast features
%        options.mask    = 5;                    % display results
%        [X2,Xn2] = Bfx_basicint(J,R,options);   
%        X3 = [X1 X2]; Xn3 =[Xn1;Xn2];
%        fprintf('\nOriginal features\n');
%        Bio_printfeatures(X3,Xn3)
%        [X4,Xn4] = Bfs_nobackground(X3,Xn3);     % delete contrast features
%        fprintf('\nSelected features\n');
%        Bio_printfeatures(X4,Xn4)
%
% D.Mery, PUC-DCC, May 2012
% http://dmery.ing.puc.cl

function [f_new,fn_new] = Bfs_nobackground(f,fn)


[ix,fn_new] = Bio_findex(fn,'contrast',0);
f_new = f(:,ix);
