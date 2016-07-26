% [f_new,fn_new] = Bfs_norotation(f,fn)
%
% Toolbox: Balu
%     This procedure deletes all no rotation invariant features.
%     It deletes the features that have in their name the strings:
%     - orient             
%     - Gabor(             
%     - [8,u2] for LBP             
%     - sLBP          
%
% Example:
%       options.b = Bfx_build({'haralick','lbp'});
%       options.colstr = 'rgb';                                 % R image
%       I = imread('testimg1.jpg');                             % input image
%       R = Bim_segbalu(I);                                     % segmentation
%       [X,Xn] = Bfx_int(I,R,options);                          % intensity features
%       [X1,Xn1] = Bfs_norotation(X,Xn);
%       Bio_printfeatures(X1,Xn1)
%
%
% See also Bfs_noposition, Bfs_nobackground.
%
% D.Mery, PUC-DCC, May 2012
% http://dmery.ing.puc.cl

function [f_new,fn_new] = Bfs_norotation(f,fn)

f_new = f; fn_new = fn;

[ix,fn_new] = Bio_findex(fn_new,'Orientation',0); f_new = f_new(:,ix);
[ix,fn_new] = Bio_findex(fn_new,'Ellipse-orient',0); f_new = f_new(:,ix);
[ix,fn_new] = Bio_findex(fn_new,'Gabor(',0); f_new = f_new(:,ix);
[ix,fn_new] = Bio_findex(fn_new,'[8,u2]',0); f_new = f_new(:,ix);
[ix,fn_new] = Bio_findex(fn_new,'sLBP'  ,0); f_new = f_new(:,ix);
[ix,fn_new] = Bio_findex(fn_new,'Fourier Abs ('  ,0); f_new = f_new(:,ix);
[ix,fn_new] = Bio_findex(fn_new,'Fourier Ang ('  ,0); f_new = f_new(:,ix);
[ix,fn_new] = Bio_findex(fn_new,'DCT('  ,0); f_new = f_new(:,ix);

