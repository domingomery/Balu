% Bex_ellipses
%
% Toolbox: Balu
%    Example 1: Ellipse estimation of regions.
%    Example 2: Selection of rices according to orientation and size.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl

clt
b(1).name = 'fitellipse';    b(1).options.show=0;  % ellipse fetaures
op.b = b;

e = input('Forms = 1, Rices = 2? ');

if e==1
    % Input Image
    I = imread('forms.bmp');

    % Segmentation
    J = imdilate(not(I),ones(3,3));
    figure(1)
    clf
    imshow(J);
    hold on
    [L,n] = bwlabel(J,4);

    % Feature extraction
    [X,Xn] = Bfx_geo(L,op);

    % Output
    Bio_drawellipse(X,'b');
else
    % Input image
    I = double(imread('rice.png'));
    figure(1)
    imshow(I,[])
    title('Original image');

    % Segmentation
    [N,M] = size(I);
    K = zeros(N,M);
    for i=1:N
        K(i,:) = I(i,:)-min(I(i,:))*ones(1,M);
    end
    T = bwareaopen(K>60,20);
    J = imclearborder(T);
    figure(2)
    imshow(J)
    title('Segmented image');
    [L,n] = bwlabel(J,4);

    % Feature extraction
    disp('extracting features...')
    [X,Xn] = Bfx_geo(L,op);

    % Output
    figure(1)
    Bio_drawellipse(X,'b');
    enterpause

    figure(3)
    k = 7;
    jj = X(:,k)<300;
    hist(X(jj,k))
    hold on
    title(['histogram for ' Xn(k,:)]);
    disp('Too small and too large rices are eliminated...');
    m1 = 150; m2 = 200;
    plot([m1 m1],[0 50],'r');
    plot([m2 m2],[0 50]),'r';
    text(m2+2,40,'big rices')
    text(m1-20,40,'small rices')
    ii = find(X(:,k)>m1&X(:,k)<m2);
    R = zeros(N,M);
    for r=1:length(ii)
        R = or(R,L==ii(r));
    end
    figure(4)
    imshow(R)
    title('selected rices')
    while(1)

        theta = input('orientation [grad]? ');
        t1 = (theta-10)*pi/180;
        t2 = (theta+10)*pi/180;
        ka = 5;
        ii = find(X(:,ka)>t1&X(:,ka)<t2&X(:,k)>m1&X(:,k)<m2);
        S = zeros(N,M);
        for r=1:length(ii)
            S = or(S,L==ii(r));
        end
        figure(5)
        imshow(S)
        title(sprintf('oriented rices at %d [grad]',theta))
    end

end

