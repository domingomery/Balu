% [f_new,fn_new] = Bfs_noposition(f,fn)
%
% Toolbox: Balu
%     This procedure deletes the features related to the position.
%     It deletes the following features:
%        - center of grav i
%        - center of grav j
%        - Ellipse-centre i
%        - Ellipse-centre j
%
% Example:
%      I = imread('testimg1.jpg');            % input image
%      R = Bim_segbalu(I);                    % segmentation
%      [X1,Xn1] = Bfx_basicgeo(R);            % basic geometric features
%      [X2,Xn2] = Bfx_fitellipse(R);          % Ellipse features
%      X3 = [X1 X2]; Xn3 =[Xn1;Xn2];
%      fprintf('\nOriginal features\n');
%      Bio_printfeatures(X3,Xn3)
%      [X4,Xn4] = Bfs_noposition(X3,Xn3);     % delete position features
%      fprintf('\nSelected features\n');
%      Bio_printfeatures(X4,Xn4)
%
%  See also Bfs_norotation, Bfs_nobackground.
%
% D.Mery, PUC-DCC, Nov. 2009
% http://dmery.ing.puc.cl

function [f_new,fn_new] = Bfs_noposition(f,fn)




s = ['center of grav i [px]   '
     'center of grav j [px]   '
     'Ellipse-centre i [px]   '
     'Ellipse-centre j [px]   '];
[n,m] = size(s);

M = size(fn,1);

f_new = f;
fn_new = fn;

ii = [];
for i=1:n;
    D = sum(abs(fn(:,1:m)-ones(M,1)*s(i,:)),2);
    ii = [ii; find(D==0)];
end
if not(isempty(ii))
    f_new(:,ii)  = [];
    fn_new(ii,:) = [];
end


