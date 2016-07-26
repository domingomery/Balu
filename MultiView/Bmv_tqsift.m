% D = Bmv_tqsift(Iq,It,options)
%
% Toolbox: Balu
%
%    Search of image query (Iq) in target image (It) using SIFT.
%
%    options.q    : sliding windows's size in pixels
%    options.d    : sliding step in pixels
%    options.nkp  : minimal number of matching keypoints
%    options.fast : '1' computes all SIFT keypoints of It at once
%                   '0' computes the SIFT keypoints for each sliding window
%    options.show : display results
%    options.roi  : region of interest where the matching in target image
%                   will be searched. If roi is not given, it will be
%                   considered that the search are is the whole image It.
%
%    D is the detection map.
%
%
%    Example 1:
%       I = imread('X1.png');
%       Iq = I(165:194,80:109);
%       It = imread('X2.png');
%       op.q=30; op.d=5; op.show=1;op.nkp=2;op.fast=1;
%       D = Bmv_tqsift(Iq,It,op);
%       figure
%       Bio_edgeview(It,bwperim(D>18))
%
%    Example 2:
%       I = imread('X1.png');
%       ix = 315:394; jx=60:139; 
%       m1 = [mean(ix) mean(jx) 1]';
%       Iq = I(ix,jx);
%       It = imread('X2.png');
%       figure(1);
%       imshow(I,[]);
%       title('Query image: blue box')
%       hold on;
%       plot(m1(2),m1(1),'rx')
%       plot([min(jx) min(jx) max(jx) max(jx) min(jx)],[max(ix) min(ix) min(ix) max(ix) max(ix)])
%       figure(2);
%       imshow(It,[]);
%       title('Target image')
%       enterpause
%       disp('Searching matching region without epiolar restriction...')
%       op.q=80; op.d=5; op.show=1;op.nkp=3;op.fast=0;
%       D1 = Bmv_tqsift(Iq,It,op);
%       figure(3)
%       Bio_edgeview(It,bwperim(D1>10))
%       title('without epiplar line')
%       enterpause
%       disp('Searching matching region with epiolar restriction...')
%       close all
%       F  = Bmv_fundamentalSIFT(I,It);        % F matrix estimation
%       ell = F*m1;                            % epipolar line
%       R = Bmv_line2img(ell,size(It));
%       op.roi = imdilate(R,ones(20,20));
%       D2 = Bmv_tqsift(Iq,It,op);
%       figure(3)
%       Bio_edgeview(It,bwperim(D2>10))
%       title('with epiplar line')
%       hold on
%       Bmv_epiplot(F,m1);
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl


function D = Bmv_tqsift(Iq,It,options)



q    = options.q;   % sliding windows's size
d    = options.d;   % sliding step
nm   = options.nkp; % minimal number of matching keypoints
show = options.show;

[N,M] = size(It);
D = zeros(N,M);


if ~isfield(options,'roi')
    Iroi = ones(N,M);
else
    Iroi = options.roi;
end

q2 = fix(q/2);

if show
    figure(1)
    clf
    imshow(Iq,[])
    title('Query image');
    hold on
    figure(2)
    clf
    imshow(It,[])
    title('Target image');
    hold on
    drawnow
end

if size(It,3)==3
    It = single(rgb2gray(It));
else
    It = single(It);
end

if size(Iq,3)==3
    Iq = single(rgb2gray(Iq)) ;
else
    Iq = single(Iq);
end

if show
    disp('Bmv_tqsift: computing SIFT descriptors...');
end

[fq, dq] = vl_sift(Iq) ;

if ~options.fast
    for i=1:d:N-q+1
        for j=1:d:M-q+1
            if Iroi(i+q2,j+q2)
                Itij = It(i:i+q-1,j:j+q-1);
                [ftij, dtij] = vl_sift(Itij);
                matches = vl_ubcmatch(dtij,dq);
                if length(matches)>nm
                    D(i:i+q-1,j:j+q-1)=D(i:i+q-1,j:j+q-1)+1;
                    if show
                        plot(j+q2,i+q2,'go');
                        drawnow
                    end
                end
            end
        end
    end
else
    [ft, dt] = vl_sift(It) ;
    jt       = ft(1,:);
    it       = ft(2,:);
    for i=1:d:N-q+1
        ii = find(and(it>=i,it<i+q));
        dti = dt(:,ii);
        jti = jt(ii);
        for j=1:d:M-q+1
            if Iroi(i+q2,j+q2)
                dtij = dti(:,and(jti>=j,jti<j+q));
                matches = vl_ubcmatch(dtij,dq);
                if length(matches)>nm
                    D(i:i+q-1,j:j+q-1)=D(i:i+q-1,j:j+q-1)+1;
                    if show
                        plot(j+q2,i+q2,'go');
                        drawnow
                    end
                end
            end
        end
    end
end