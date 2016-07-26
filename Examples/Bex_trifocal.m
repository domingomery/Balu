clt
f.path          = '/Volumes/domingomery/Mingo/Matlab/balu3/';
f.prefix        =  'X';   f.extension   =  '.png';
f.digits        = 1;
f.gray          = 1;
f.subsample     = 1;
f.resize        = 1;
f.imgmin        = 1;
f.imgmax        = 6;
f.window        = [];
f.negative      = 0;
f.sequence      = 1:f.imgmax;

show = 0;                        % use show = 1 to display results
op1.show        = show;          % 1: display results
op1.descriptor  = 'harris+sift'; % keypoints' descriptors
op1.sfm_method  = 2;             % SfM > 1: affine, 2: projective (see Btr_sfm)
op1.sfm_samples = 4;             % all samples for RANSAC SfM (see Btr_sfm)
op1.sfm_iter    = 100;           % iterations for RANSAC SfM (see Btr_sfm)
op1.sfm_cover   = 0.5;           % cover of the SfM solution (see Btr_sfm)
op1.matching    = 1;             % 1: direct, 2: homomography (see Btr_sift2)
op1.sort        = 1;             % sort image sequence using Bsq_sort
op1.indexonly   = 0;             % sort and change image sequence (see Bsq_sort)

[P,fn,kp1] = Btr_sfseq(f,op1);
figure(4);Bsq_show(fn,6);

while (1)
    i = input('image 1? ');
    j = input('image 2? ');
    k = input('image 3? ');

    A = P(indices(i,3),:);
    B = P(indices(j,3),:);
    C = P(indices(k,3),:);

    T = Bmv_trifocal(A,B,C);
    F12 = Bmv_fundamental(A,B);
    F13 = Bmv_fundamental(A,C);
    F23 = Bmv_fundamental(B,C);

    I1 = Bio_loadimg(fn,i);figure(1);clf;imshow(I1,[]);title('image 1');hold on; axis on
    I2 = Bio_loadimg(fn,j);figure(2);clf;imshow(I2,[]);title('image 2');hold on; axis on
    I3 = Bio_loadimg(fn,k);figure(3);clf;imshow(I3,[]);title('image 3');hold on; axis on

    new = 1;
    while(new)
        figure(1);
        disp('Select point in image 1...');
        p1 = vl_click; m1 = [p1([2 1]); 1];
        plot(m1(2),m1(1),'rx');
        figure(2);
        Bmv_epiplot(F12,m1);
        disp('Select corresponding point in image 2...');
        p2 = vl_click; m2 = [p2([2 1]); 1];
        m3s = Bmv_reproj3(m1,m2,T,2)
        plot(m2(2),m2(1),'rx');
        figure(3)
        Bmv_epiplot(F13,m1);
        Bmv_epiplot(F23,m2);
        plot(m3s(2),m3s(1),'w*');
        plot(m3s(2),m3s(1),'wo');
        new = input('New points?');
    end
end









