%  Dt = Bim_regiongrow(I,th)
%
% Toolbox: Balu
%  Interactive selection of regions with similar grayvalues using growing
%  region algorithm.
%  I   : grayscale image
%  th  : maximal difference of grayvalues in the region
%        (default: th = 20)
%  Dt  : output binary image with segmented regions
%
%  Example:
%     I = imread('X1.png');
%     Dt = Bim_regiongrow(I);
%
% D.Mery, PUC-DCC, Apr. 2011
% http://dmery.ing.puc.cl
%

function Dt = Bim_regiongrow(I,th)
I = double(I);
Dt = zeros(size(I));
op = 1;
points=1;
if ~exist('th','var')
    th = (max2(I)-min2(I))/16;
end
while(points)
    D = zeros(size(I));
    N = ones(size(I));
    if op==1
        figure(2)
        imshow(Dt)
        title('Selected regions');
        figure(1)
        Bio_edgeview(I,bwperim(Dt));
        disp('select pixel of a new region in Figure 1... (or press Esc to end the selection)')
        p = vl_click;
        if isempty(p)
            break
        end
        i = round(p(2));
        j = round(p(1));
        m = I(i,j);

    end
    D(i,j) = 1;
    ok = 0;
    nd0 = sum2(D);
    while not(ok)
        nd = nd0;
        %jj = find(D==1);
        %m = mean(I(jj));
        Dd = and(N,imdilate(D,ones(3,3)));
        S  = Dd-D;
        ii = find(S==1);
        for k=1:length(ii)
            if abs(I(ii(k))-m)<th
                D(ii(k))=1;
            else
                N(ii(k))=0;
            end
        end
        nd0 = sum2(D);
        ok = nd0-nd==0;
    end
    figure(1)
    Bio_edgeview(I,bwperim(or(Dt,D)));
    figure(2)
    imshow(D)
    disp(' ')
    disp('1: accept this region and take another');
    disp('2: delete this region and take another');
    disp('3: the same region but smaller');
    disp('4: the same region but larger');
    disp('5: cut this region');
    disp('6: accept this and end, no more regions');
    disp('7: delete this and end, no more regions');
    disp(' ')
    op = input('Select option: ');
    switch op
        case 1
            Dt = or(D,Dt);
        case 2
            op=1;
        case 3
            th = max(th-2,0);
        case 4
            th = th+2;
        case 5
            ok = 0;
            while not(ok)
                figure(2)
                imshow(D)
                figure(1)
                clf
                Bio_edgeview(I,bwperim(or(Dt,D)));
                hold on
                disp('Define Bounding box to cut in Figure 1 with two clicks...')
                R = ones(size(I));
                p = round(vl_click(2));
                p_i = sort(p(2,:));
                p_j = sort(p(1,:));
                plot(p_j([1 1 2 2 1]),p_i([1 2 2 1 1]))
                R(p_i(1):p_i(2),p_j(1):p_j(2)) = 0;
                Do = and(R,D);
                figure(2)
                imshow(Do)
                op = input('1 = ok, 2 = delete this region, 3 = extra cut, 4 = repeat this cut, 5 = accept whole box? ');
                switch op
                    case 1
                        ok = 1;
                    case 2
                        ok = 1;
                        Do = zeros(size(I));
                    case 3
                        D = Do;
                    case 5
                        Do(p_i(1):p_i(2),p_j(1):p_j(2))=1;
                        ok = 1;
                end
            end
            D = Do;
            Dt = or(D,Dt);
            op = 1;
            figure(1)
            clf
        case 6
            points=0;
            Dt = or(D,Dt);
        case 7
            points=0;
    end
end
figure(2)
imshow(Dt)
figure(1)
Bio_edgeview(I,bwperim(Dt));

