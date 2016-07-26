% Bprintfeatures(X,Xn)    % for features with feature names
% Bprintfeatures(X)       % for feature values only
%
% Toolbox: Balu
%    Display extracted features.
%    Xn: feature names  (matrix mxp, one row per each string of p
%        characters)
%    X:  feature values (vector 1xm for m features)
%    Xu: feature units  (matrix mxq, one row per each string of q
%        characters)
%   
%    These variables are the outputs of Bgeofeatures or Bintfeatures.
%    
%    The output of Bprintfeatures is like this:
%
%    1 center of grav i       [pixels]         163.297106
%    2 center of grav j       [pixels]         179.841850
%    3 Height                 [pixels]         194.000000
%    4 Width                  [pixels]         196.000000
%    5 Area                   [pixels]         29361.375000
%    :         :           :                  :
%
%   Example 1: Display of standard geometric features of testimg1.jpg
%      I = imread('testimg1.jpg');            % input image
%      R = Bsegbalu(I);                       % segmentation
%      [X,Xn,Xu] = Bfg_standard(R);           % standard geometric features
%      Bprintfeatures(X,Xn,Xu)                 
%
%   Example 2: Display of first 5 samples of datagauss.mat
%      load datagauss
%      Xn = ['[length]';'[weigh] '];
%      Xu = ['cm';'kg'];
%      for i=1:5
%         fprintf('Sample %d:\n',i);
%         Bprintfeatures(X(i,:),Xn,Xu)
%         Benterpause
%      end
%
%
%   See also Bplotfeatures.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function Bprintfeatures(X,Xn)
N = length(X);
if ~exist('Xn','var')
    Xn = char(zeros(N,16));
end
for k=1:size(Xn,1)
    fprintf('%3d %s %f\n',k,Xn(k,:),X(k));
end    
