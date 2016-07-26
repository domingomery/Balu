% [X,Xn] = Bfx_geo(R,options)
% [X,Xn] = Bfx_geo(L,options)
%
% Toolbox: Balu
%
%    Gemteric feature extraction.
%
%    This function calls gemetric feature extraction procedures of binary
%    image R or labelled image L.
%
%    X is the feature matrix (one feature per column, one sample per row),
%    Xn is the list with the names of these features (see Example
%    to see how it works).
%
%   Example 1: Extraction of one region image
%      b(1).name = 'hugeo';       b(1).options.show=1;         % Hu moments
%      b(2).name = 'flusser';     b(2).options.show=1;         % Flusser moments
%      b(3).name = 'fourierdes';  b(3).options.show=1;         % Fourier
%                                 b(3).options.Nfourierdes=12; % descriptors
%      options.b = b;
%      I = imread('testimg1.jpg');                             % input image
%      R = Bim_segbalu(I);                                     % segmentation
%      [X,Xn] = Bfx_geo(R,options);                            % geometric features
%      Bio_printfeatures(X,Xn)
%
%   Example 2: Extraction of multiple regions image
%      b(1).name = 'hugeo';       b(1).options.show=1;         % Hu moments
%      b(2).name = 'basicgeo';    b(2).options.show=1;         % basic geometric fetaures
%      b(3).name = 'fourierdes';  b(3).options.show=1;         % Fourier
%                                 b(3).options.Nfourierdes=12; % descriptors
%      options.b = b;
%      I = imread('rice.png');                                 % input image
%      [R,m] = Bim_segmowgli(I,ones(size(I)),40,1.5);          % segmentation
%      [X,Xn] = Bfx_geo(R,options);                            % geometric features
%      figure; hist(X(:,12));xlabel([Xn(12,:)])                % area histogramm
%      ii = find(abs(X(:,20))<15);                             % rice orientation
%      K = zeros(size(R));                                     % between -15 and 15 grad
%      for i=1:length(ii);K=or(K,R==ii(i));end
%      figure; imshow(K);title('abs(orientation)<15 grad')
%
%   See also Bfx_basicgeo, Bfx_hugeo, Bfx_flusser, Bfx_gupta,
%            Bfx_fitellipse, Bfx_fourierdes, Bfx_files.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [X,Xn] = Bfx_geo(R,options)


if isempty(R)
    error('Bfx_geo: R is empty. Geometric features without segmentation has no sense.');
else
    b = options.b;
    n = length(b);
    m = int16(max(R(:)));
    X = [];
    for j=1:m
        Rj = R==j;
        Xj = [];
        Xnj = [];
        for i=1:n
            s = b(i).name;
            if sum(s(1:3)=='Bfx')~=3
                s = ['Bfx_' s];
            end
            if  ~exist(s,'file')
                error(sprintf('Bfx_geo: function %s does not exist.',b(i).name))
            end

            [Xi,Xni] = feval(s,Rj,b(i).options);
            Xj = [Xj Xi];
            Xnj = [Xnj; Xni];
        end
        X = [X;Xj];
    end
    Xn = Xnj;
end
