% From Spyrou et al 2010:
% ... most of the common visual words, i.e., those with smallest iDF values,
% are not discriminative and their absence would facilitate the retrieval
% process. On the other hand, the rarest visual words are in most of the
% cases as result of noise and may distract the retrieval process. To
% overcome those problems, a stop list is created that includes the most
% and the least frequent visual words of the image collection.
%
% From Sivic & Zisserman (2003):
% Using a stop list the most frequent visual words that occur in almost all
% images are suppressed. The top 5% (of the frequency of visual words over 
% all the keyframes) and bottom 10% are stopped.

function [Voc,ix,i_stop,i_go,kd,kd_go] = Bsq_stoplist(Y,options)

top_bound    = options.top_bound;        % eg. 0.05  top    : most frequent word to be stopped
bottom_bound = options.bottom_bound;     % eg. 0.10  bottom : less frequent word to be stopped
NV           = options.NV;               % number of visual words
ixn          = options.ixn;              % indices of Y rows per document

ii_top       = round(NV*(1-top_bound))+1:NV;
ii_bottom    = 1:round(NV*bottom_bound);
ii_stop      = [ii_bottom ii_top]; 
show         = options.show;

if isfield(options,'newdictionary');
    newdictionary = options.newdictionary;
else
    newdictionary = 1;
end

if newdictionary == 1
    fprintf('Computing vocabulary of Y (%dx%d) with %d visual words...',size(Y,1),size(Y,2),NV);
    tic
    [Voc,ix] = vl_kmeans(Y',NV,'Algorithm','Elkan');
    t = toc;
    fprintf('in %f sec.\n',t);
    
else
    fprintf('Using pre-computed vocabulary with %d visual words...\n',NV);
    Voc      = options.voc;
    ix       = options.ix_voc;
end
ww           = double(ix);
N            = max(ixn);
TF           = zeros(N,NV);

for i=1:N
    ii = ixn==i;
    ww_i = ww(ii);
    TF(i,:) = hist(ww_i,1:NV);
end
% term frequency TF(d,t) = number of occurrences of word t in document d, d=1,...k, t=1,...NV

DF = sum(TF>0); % document frequency (DF(t) = number of documents in the collection that contain a word t, t=1,...NV)

% Not necessary to create a stop list
%    iDF = log(N./(DF));
%    TFIDF = TF.*(ones(k,1)*iDF);
%    [ii,jj] = sort(iDF,'descend');

[~,jj] = sort(DF);
i_stop = jj(ii_stop);
i_go   = jj; i_go(ii_stop) = [];

kd     = vl_kdtreebuild(Voc);
kd_go  = vl_kdtreebuild(i_go);


% V_stop = Voc(:,j_stop);
% V_go   = Voc(:,j_go);

% like Fig. 5 of Sivic & Zisserman (2003)
if show
    figure(12);clf;
    subplot(1,3,1); mesh(TF);
    subplot(1,3,2);plot(sort(DF,'descend'));ax=axis;
    subplot(1,3,3);plot(sort(DF(i_go),'descend'));axis(ax)
end



    