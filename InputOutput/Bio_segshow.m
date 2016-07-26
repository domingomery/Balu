% Bio_segshow(I,p)
%
% Toolbox: Balu
%    Display original image and segmented image
%
%    Bimshow display 4 images:
%       1: Original image (matrix I)
%       2: Segmented (using command Bsegbalu)
%       3: High contrast (using command Bsegbalu)
%       4: Edges (using Bedgeview)
%
%       p is the parameter used by command Bsegbalu.
%
% Example 1: Segmentation using Balu segmentation algorithm
%     I = imread('testimg1.jpg');
%     Bio_segshow(I)
%     Repeat this examples for images testimg2, testimg3 and testimg4. Last
%     test image requires R = Bimshow(I,-0.1) for better results.
%
% Example 2: Segmentation using PCA and with more sensibility
%     I = imread('testimg2.jpg');
%     Bio_segshow(I,'Bim_segpca',-0.15)
%
% Example 3: The name of the image can be given as argument
%     Bio_segshow('testimg4.jpg','Bim_segmaxfisher',-0.1)
%
% See also Bim_segbalu.
%
% D.Mery, PUC-DCC, May 2008
% http://dmery.ing.puc.cl
%
function Bio_segshow(I,segname,p)


if not(exist('segname','var'))
    segname = 'Bim_segbalu';
end


if compare(class(I),'char')==0
    I = imread(I);
end


if not(exist('p','var'))
    p = 0;
end
fsname = ['[R,E,J] = ' segname '(I,p);'];
eval(fsname);


subplot(2,2,1);imshow(I);title('original image')
subplot(2,2,2);imshow(R);title('segmented image')
subplot(2,2,3);imshow(J);title('high contrast image')
subplot(2,2,4);Bio_edgeview(I,imdilate(E,ones(3,3)));title('edge image')
