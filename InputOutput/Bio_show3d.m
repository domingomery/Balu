function Bio_show3d(J,images,map,sl)
[N,M,P,K]=size(J);
if K>1
    H = zeros(N,M,K);
    for i=1:K
        H(:,:,i) = J(:,:,1,i);
    end
    J = H;
    P = K;
end

if ~exist('images','var')
    images = [];
end

if ~exist('sl','var')
    sl = 3;
end


if ~exist('map','var')
    map = [];
end
%J = uint8(round(J));

switch sl
    case 1
        if isempty(images)
            images = 1:N;
        end
        G = zeros(M,P);
        for i=images;
            G(:) = J(i,:,:);
            imshow(G',map);
            drawnow
            enterpause(0)
        end
    case 2
        if isempty(images)
            images = 1:M;
        end
        G = zeros(N,P);
        for i=images;
            G(:) = J(:,i,:);
            imshow(G',map);
            drawnow
            enterpause(0)
        end
    case 3
        if isempty(images)
            images = 1:P;
        end
        for i=images;
            imshow(J(:,:,i),map);
            drawnow
            enterpause(0)
        end
end