% [v,Xcen,H] = Bsq_vocabulary(kp,V,options)
%
% Toolbox: Balu
%
%    Visual vocabulary.
%
%    kp.des is the descriptions of the keypoint.
%    kp.img is a vector containing the number of the image of each keypoint.
%    The minimum value of img is allways 1 and the maximum is the number of
%    processed images, i.e., f.imgmax-f.imgmin+1.
%    V is the number of visual words (clusters).
%    options.show display results.
%
%    v is the document representation using "term frequency-inverse
%    document frequency"
%    Xcen are the centroids of the kmeans.
%    H are the histograms.
%
%
%    Reference:
%      Sivic & Zisserman: Efficient visual search for videos cast as text
%      retrieval. 31(4):591-606. PAMI 2009.
%
%    Example:
%    f.path             = ''; % Balu directory as path or current directory
%    f.extension        = '.jpg';
%    f.prefix           = 'testimg';
%    f.digits           = 1;
%    f.gray             = 1;
%    f.subsample        = 1;
%    f.resize           = [256 256];
%    f.window           = [];
%    f.negative         = 0;
%    f.sequence         = 1:6;
%    f.imgmin           = 1;
%    f.imgmax           = 6;
%    options.show       = 1;
%    options.descriptor = 'sift';
%    options.clean      = 0;
%    kp = Bsq_des(f,options);
%    [v,Xcen,H] = Bsq_vocabulary(kp,100,options);
%
% See also Bsq_vgoogle, Bsq_sort.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [v,Xcen,H,Vnew,kd,N,Ni,Xcen_old,jsel] = Bsq_vocabulary(kp,V,options)

if ~exist('options','var')
    options.show = 0;
end

if ~isfield(options,'tfidf')
    options.tfidf = 1;
end

if ~isfield(kp,'img')
    kp.img = ones(size(kp.des),1);
end

if ~isfield(options,'show')
    options.show = 0;
end

show = options.show;

if show
    fprintf('Bsq_vocabul: Vector Quantization with %d clusters...\n',V)
end
X  = single(kp.des);
Xcen = (vl_kmeans(X',V,'Algorithm','ANN'))';
Xcen_old = Xcen;
if show
    disp('Bsq_vocabul: building kd-tree...')
    pause(0)
end

kd = vl_kdtreebuild(Xcen');

H    = [];
Vnew = V;
v    = [];


if options.tfidf == 1



    N = max(kp.img);

    H = zeros(N,V);
    Nd = zeros(N,1);
    for d=1:N
        h = zeros(1,V);
        if show
            fprintf('Bsq_vocabul: processing image %d\n',d);
        end
        ii = find(kp.img==d);
        if ~isempty(ii)
            Xt = single(kp.des(ii,:));
            j  = vl_kdtreequery(kd,Xcen',Xt','NumNeighbors',1)';
            nd = length(j);
            Nd(d) = nd;
            for k=1:nd;
                h(j(k))=h(j(k))+1;
            end
            H(d,:) = h;
        end
    end


    H1 = H>0;
    hs = mean(H1);
    j = hs>0.85;
    H(:,j) = [];
    jsel = not(j);
    Xcen(j,:) = [];
    V = size(H,2);

    v = zeros(N,V);

    Ni  = sum(H>0,1);
    for d=1:N
        nd  = Nd(d);
        if (nd>0)
            for i=1:V
                nid = H(d,i);
                ti = nid/nd*log(N/Ni(i));
                v(d,i) = ti;
            end
        end
    end

    Vnew = V;
end