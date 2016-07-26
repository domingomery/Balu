% From Spyrou et al 2010:
% ... most of the cmmon visual words, i.e., those with smallest iDF values,
% are not descriminative and their abscence would facilitate the retrieval
% process. On the other hand, the rarest visual words are in most of the
% cases as result of noise and may distract the retrieval process. To
% overcome those problems, a stop list is created that includes the most
% and the least frequent visual words of the image collection.
%
% From Sivic & Zisserman (2003):
% Using a stop list the most frequent visual words that occur in almost all
% images are supressed. The top 5% (of the frequency of visual words over 
% all the keyframes) and bottom 10% are stopped.
%
% Term Frequency-Inverse Document Frequency Weighting
% After Sivic & Zisserman (PAMI 2009). See Section 4.1
% Each image of the training database is one 'document'
% The words of each document will be filtered out using a stop list
% (see program asr_stoplist.m).
%
% D.Mery, Notre Dame (2014)

function opvoc = Bsq_visualvoc(Y,options)

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
    if show>0
    fprintf('     computing visual vocabulary from (%dx%d) with %d visual words...',size(Y,1),size(Y,2),NV);
    tic
    end
    % [opvoc.voc,opvoc.ix_voc] = vl_kmeans(Y',NV,'Algorithm','Elkan');
    [opvoc.voc,opvoc.ix_voc] = vl_kmeans(Y',NV,'Algorithm','Elkan','NumRepetitions',3);
    if show>0
    t = toc;
    fprintf('in %f sec.\n',t);
    end
else
    fprintf('Using pre-computed vocabulary with %d visual words...\n',NV);   
end

ww           = double(opvoc.ix_voc);
N            = max(ixn);                   % # documents in the whole database
nid          = zeros(N,NV);                % nid(d,i) # occurences of word i in document d 

ii           = 1:NV;                       % indices of the words of the vocabulary

for d=1:N
    dd = ixn==d;                           % indices of words in Y of document d
    ww_d = ww(dd);
    nid(d,:) = hist(ww_d,ii);              % # occurrences of word i=1...NV in document d
end
% term frequency nid(d,i) = number of occurences of word i in document d, d=1,...k, t=1,...NV

Ni = sum(nid>0); % document frequency (Ni(i) = # of documents in the collection containing word i, i=1,...NV)


[sNi,jj]              = sort(Ni);
opvoc.i_stop        = jj(ii_stop);
opvoc.i_go          = jj; 
opvoc.i_go(ii_stop) = [];

opvoc.i_stop        = sort(opvoc.i_stop);
opvoc.i_go          = sort(opvoc.i_go);


opvoc.kd_voc           = vl_kdtreebuild(opvoc.voc);
opvoc.kd_go            = vl_kdtreebuild(opvoc.i_go);
opvoc.voc_stop         = opvoc.voc(:,opvoc.i_stop);
opvoc.voc_go           = opvoc.voc(:,opvoc.i_go);
opvoc.Ni               = Ni;
Ni(opvoc.i_stop)       = 0;
opvoc.Ni_go            = Ni;

nid_go                   = nid;
nid_go(:,opvoc.i_stop) = 0;

if options.tfidf > 0
% The following steps are necessary to create a stop list
% the frequencies are computed after eliminating the stop words
    nd            = sum(nid_go,2); % nd(d): # documents in document d (stop list words are not included)
    iDF           = log(double(N)./(double(Ni)+1e-20));
    opvoc.Vd    = Bft_uninorm((nid_go./(nd*ones(1,NV))).*(ones(N,1)*iDF));
end

% similar to Fig. 5 of Sivic & Zisserman (2003)
if show>2
    figure(12);clf;
    subplot(1,2,1); mesh(nid);
    subplot(1,2,2);plot(sNi);ax=axis;
    hold on
    plot([min(ii_top)    min(ii_top)]   ,[ax(3) ax(4)],'r:')
    plot([max(ii_bottom) max(ii_bottom)],[ax(3) ax(4)],'r:')
end



    