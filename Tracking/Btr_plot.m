% Btr_plot(kp,A,files,options)
%
% Toolbox: Balu
%
%    Plot of tracks.
%
%    kp keypoints structure according function Bsq_des (see help)
%
%    A is a Nxm matrix with N matchings in m views. Each row is a track to be plot.
%
%    files is a structure that define the images of the sequence according
%    to function Bio_loadimg (see help).
%
%    options.plottraj   : 1 plot the trajectory
%    options.plotimg    : 1 display the image number
%    options.plotpoints : 1 display the index number of the keypoint 
%    options.img1       : 1 first image to be displayed
%    options.plotsquare : 1 plot a bounding box arroind the keypoint 
%
%  Example:
%      See example in Btr_demo.
%
%  See also Bsq_des, Bio_loadimg, Btr_demo.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function Btr_plot(kp,A,files,options)

f          = files;
img        = kp.img;
fra        = kp.fra;
ilu        = kp.ilu;
traj       = options.plottraj;
imgshow    = options.plotimg;
pointshow  = options.plotpoints;
i1         = options.img1;
squareshow = options.plotsquare;


if isempty(A)
    A = [(1:length(img))' (1:length(img))'];
end


if i1 == 0;
    n = size(A,1);
    if n==1;
        options.img1 = img(A(1,1));
        Btr_kplot(kp,A,files,options)
    else
        if n>1
            p1 = 1;
            p2 = f.imgmax-f.imgmin+1;
            for p=p1:p2
                options.img1 = p;
                Btr_plot(kp,A,files,options)
            end
        end

    end
else

    ii = find((img(A(:,1))==i1));
    n = length(ii);
    if n>0
        if traj==0
            for i=1:n
                Y = A(ii(i),:);
                Y = Y(Y>0);
                x = fra(Y,1:3);
                nx = size(x,1);
                figure(1)
                clf
                for j=1:nx
                    imgj = ilu(img(Y(j)));
                    X = Bio_loadimg(f,imgj);
                    subplot(1,nx,j)
                    imshow(X,[])
                    hold on
                    title(sprintf('image %d',imgj))
                    if imgshow
                        text(x(j,1)-5,x(j,2)+9,num2str(imgj))
                    end
                    if pointshow
                        text(x(j,1)-5,x(j,2)-9,num2str(Y(j)))
                    end
                    vl_plotframe([x(j,1:2) x(j,3)*3]);
                    plot(x(j,1),x(j,2),'rx')
                    pause(0.5)

                end
            end
        else
            X = Bio_loadimg(f,ilu(i1));
            clf
            imshow(X,[])
            hold on
            title(sprintf('Image %d',ilu(i1)))
            for i=1:n
                Y = A(ii(i),:);
                Y = Y(Y>0);
                x = fra(Y,1:2);
                s = 2*fra(Y,3);
                plot(x(:,1),x(:,2))
                plot(x(:,1),x(:,2),'o')
                for j=1:size(x,1)
                    if imgshow
                        text(x(j,1)-5,x(j,2)+9,num2str(img(Y(j))))
                    end
                    if pointshow
                        text(x(j,1)-5,x(j,2)-9,num2str(Y(j)))
                    end
                    if squareshow
                        x1 = x(j,1);
                        y1 = x(j,2);
                        a  = s(j);
                        xx = [x1-a x1-a x1+a x1+a x1-a];
                        yy = [y1-a y1+a y1+a y1-a y1-a];
                        plot(xx,yy,'r')
                    end
                end
            end
        end
    else

        X = Bio_loadimg(f,ilu(i1));
        clf
        imshow(X,[])
        title(sprintf('Image %d',i1))

    end


end
pause(0.1)
drawnow