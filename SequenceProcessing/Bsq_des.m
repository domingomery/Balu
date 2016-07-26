% kp = Bsq_des(f,options)
%
% Toolbox: Balu
%
%    Description of a sequence.
%
%    f is the structure that defines the sequence
%
%    f.path          : directory where are the files
%    f.extension     : extension (eg: 'jpg')
%    f.prefix        : prefix (eg: 'DSC_')
%    f.digits        : number of digits (eg:4)
%    f.gray          : 1 means rgb to gray conversion
%    f.subsample     : subsampling rate (eg:1 means no subsample)
%    f.resize        : parameter of imresize, 0 means no imresize
%    f.window        : image window
%    f.negative      : negative window
%    f.sequence      : if seq = [3 2 1], the image for i=1 will be No. 3
%    f.imgmin        : minimal number of the images (not used by Bio_loadimg)
%    f.imgmax        : maximal number of the images (not used by Bio_loadimg)
%
%    options.descriptor can be: 'sift', 'surf', 'sift-plus', 'phow',
%                'harris', 'harris-defects', 'harris+sift', 'mser', 'clp',
%                'blops', 'clp+harris'.
%
%    options.show display results
%    options.param are paremeters of the method (if any).
%    options.clean is the name of the function used to eliminate keypoints
%    outside of the object of interest, this function is a tailored
%    function and returns a binary image with '1' object and '0' background
%    options.clean = 0 means no cleaning
%
%    kp (keypoints) is a structure with the following fields:
%    kp.fra and kp.des are the frames and descriptions of each keypoint.
%    kp.img is a vector containing the number of the image of each keypoint.
%    The minimum value of img is allways 1 and the maximum is the number of
%    processed images, i.e., f.imgmax-f.imgmin+1
%    kp.ilu is a look up table of the real number of the images (see example).
%
%    Example:
%    f.path             = ''; % Balu directory as path or current directory
%    f.prefix           = 'testimg';
%    f.extension        = '.jpg';
%    f.digits           = 1;
%    f.gray             = 1;
%    f.subsample        = 1;
%    f.resize           = [256 256];
%    f.window           = [];
%    f.negative         = 0;
%    f.sequence         = 1:6;
%    f.imgmin           = 5;
%    f.imgmax           = 6;
%    options.show       = 1;
%    options.descriptor = 'sift';
%    options.clean      = 0;
%
%    kp = Bsq_des(f,options);
%
%    % The frames and description of image 5 are:
%
%    i = find(kp.ilu==5);
%    j = find(kp.img==i);
%    f5 = kp.fra(j,:);
%    d5 = kp.des(j,:);
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function kp = Bsq_des(f,options)


imin  = f.imgmin;        % first image
imax  = f.imgmax;        % last image

i0 = 0;

img = [];
fra = [];
des = [];
ilu = [];


show   = options.show;
method = lower(options.descriptor);


for i=imin:imax
    if show
        fprintf('Bsq_des    : processing image %5d...',i)
    end
    J = Bio_loadimg(f,i);
    if length(J(:))>1
        i0 = i0+1;
        im = single(J);
        switch method
            case 'sift'
                [frames, descrs] = vl_sift(im) ;

            case 'surf'
                Options.verbose = false;
                Options.tresh=0.5;
                Options.octaves = 5;
                Options.init_sample = 2;
                Ipts=OpenSurf(J,Options);
                x = cat(1,Ipts.x);
                y = cat(1,Ipts.y);
                s = cat(1,Ipts.scale);
                o = cat(1,Ipts.orientation);
                frames = [x';y';s';o'];
                n = length(Ipts);
                descrs = zeros(64,n);
                descrs(:) = cat(1,Ipts.descriptor);
                % l = cat(1,Ipts.laplacian);
                % ii = l==0;
                ii = frames(3,:)<4;
                frames = frames(:,ii);
                descrs = descrs(:,ii);

            case 'sift-plus'
                H =  vl_harris( vl_imsmooth(J,3),7);
                idx = vl_localmax(H,0.3) ;
                [y,x] = ind2sub(size(im),idx);
                n = length(x);
                frames = [x; y; 20*ones(1,n);zeros(1,n)];
                % [fhar2, dhar2] = vl_sift(im,'frames',frames) ;
                % [fsift, dsift] = vl_sift(im,'EdgeThresh',50) ;
                % frames = [fhar2 fsift];
                % descrs = [dhar2 dsift];
                [fhar2, dhar2] = vl_sift(im,'frames',frames) ;
                %[fsift, dsift] = vl_sift(im,'EdgeThresh',50) ;
                frames = [fhar2];
                descrs = [dhar2];

            case 'sifts'
                [frames, descrs] = BsiftmatchS(J,options.param.Iref,options.param.t1,options.param.t2,0) ;
            case 'phow'
                [frames, descrs] = vl_phow(im);

            case 'harris'
                H =  vl_harris( vl_imsmooth(J,3),3);
                idx = vl_localmax(H ) ;
                [y,x] = ind2sub( size(im), idx );
                n = length(x);
                frames = [x; y; 20*ones(1,n) ; zeros(1,n)];
                [frames, descrs] = vl_sift(im,'frames',frames) ;

            case 'harris-defects'
                % p       = [3 0.2 5 0.5 3];
                p = [13 2 5 2 13];
                J1       = vl_imsmooth(J,p(5));
                im      = single(J1);
                H       = vl_harris(J1,p(1));
                idx     = vl_localmax(H,p(2));
                [y1,x1] = ind2sub(size(im),idx);
                H       = vl_harris(J,p(3));
                idx     = vl_localmax(H,p(4));
                [y2,x2] = ind2sub(size(im),idx);
                x       = [x1 x2];
                y       = [y1 y2];
                n       = length(x);
                frames = [x; y; 5*ones(1,n) ; zeros(1,n)];
                [frames, descrs] = vl_sift(im,'frames',frames) ;

            case {'harris+sift','sift+harris'}
                H =  vl_harris( vl_imsmooth(J,3),3);
                idx = vl_localmax(H) ;
                [y,x] = ind2sub(size(im),idx);
                n = length(x);
                frames = [x; y; 20*ones(1,n);zeros(1,n)];
                [fhar2, dhar2] = vl_sift(im,'frames',frames);
                [fsift, dsift] = vl_sift(im) ;
                frames = [fhar2 fsift];
                descrs = [dhar2 dsift];

            case 'mser'
                [r,frames] = vl_mser(uint8(J),'MinDiversity',0.7,'MaxVariation',0.2,'Delta',10) ;
                % [r,frames] = vl_mser(uint8(J),'MinDiversity',options.param(1),'MaxVariation',options.param(2),'Delta',options.param(3)) ;
                frames = vl_ertr(frames) ;
                [frames, descrs] = vl_sift(im,'frames',frames(1:4,:));

            case 'segmser'
                frames = Bim_segmser(J,options.param);
                [frames, descrs] = vl_sift(im,'frames',frames);

            case 'razor'
                frames = RazorDetector(J,options.param);
                if ~isempty(frames)
                    [frames, descrs] = vl_sift(im,'frames',frames);
                end
            case 'gun'
                frames = GunDetector(J,options.param);
                if ~isempty(frames)
                    [frames, descrs] = vl_sift(im,'frames',frames);
                end

            case 'clp'
                [Y,feat] = Bim_segdefects(J);
                n = size(feat,1);
                fclp = [feat(:,[2 1 3]) zeros(n,1)]';
                [fclp, dclp] = vl_sift(im,'frames',fclp);
                frames = fclp;
                descrs = dclp;

            case 'blops'
                % options.param = [13 20 0 50 15 200 0.5 1.2]; % para la gillette
                % options.param = [17 15 0 100 15 2000 0 0]; % para la manzana
                [Y,feat] = Bim_segblops(J,options.param);
                n = size(feat,1);
                fblops = [feat(:,[2 1 3]) zeros(n,1)]';
                [fblops, dblops] = vl_sift(im,'frames',fblops);
                frames = fblops;
                descrs = dblops;


            case 'clp+harris'
                H =  vl_harris( vl_imsmooth(J,3),3);
                idx = vl_localmax(H) ;
                [y,x] = ind2sub(size(im),idx);
                n = length(x);
                frames = [x; y; 10*ones(1,n);zeros(1,n)];
                [fhar1, dhar1] = vl_sift(im,'frames',frames) ;
                [Y,feat] = segdefects(J);
                n = size(feat,1);
                fclp = [feat(:,[2 1 3]) zeros(n,1)]';
                [fclp, dclp] = vl_sift(im,'frames',fclp) ;
                [fsift, dsift] = vl_sift(im) ;
                frames = [fclp fhar1 fsift];
                descrs = [dclp dhar1 dsift];
        end

        if ~isfield(options,'clean')
            options.clean = 0;
        end


        if sum(options.clean) > 0
            R = feval(options.clean,J);
            s = size(J);
            x = round(frames(1,:))';
            y = round(frames(2,:))';
            ix = sub2ind(s,y,x);
            t  = R(ix);
            j = t==1;
            frames = frames(:,j);
            descrs = descrs(:,j);
        end
        if show
            clf
            imshow(J,[]);
            pause(0);
            x = frames(1,:);
            y = frames(2,:);
            hold on;
            plot(x,y,'.')
            fprintf('... %6d keypoints\n',length(x(:)))
            drawnow
        end
        if ~isempty(frames)
            fra = [fra; frames'];
            des = [des; descrs'];
            img = [img; i0*ones(size(frames,2),1)];
            ilu = [ilu; i];
        end

    else
        error('Bsq_des    : image not found.');
    end
end

kp.fra = fra;
kp.des = des;
kp.img = img;
kp.ilu = ilu;

