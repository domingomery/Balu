% J = Bim_colorconv(X,s1,s2)
%
% Toolbox: Balu
%    Color conversion of image (in color space s1) into image Y
%   (in color space s2).
%
%
%  Example:
%     I = imread('testimg2.jpg');
%     J = Bim_colorconv(I,'rgb','hcm');
%     figure(1)
%     imshow(I); title('control image')
%     figure(2)
%     imshow(J); title('high contrast image')
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl


function Y = Bim_colorconv(X,s1,s2)

s1 = lower(s1);
s2 = lower(s2);

switch s1
    case 'rgb'
        switch s2
            case 'rgb'
                Y = X;
            case 'xyz'
                Y = vl_rgb2xyz(X);
            case 'lab'
                Z = vl_rgb2xyz(X);
                Y = vl_xyz2lab(Z);
            case 'luv'
                Z = vl_rgb2xyz(X);
                Y = vl_xyz2luv(Z);
            case 'hsv'
                Y = rgb2hsv(X);
            case 'hsi'
                Y = Bim_rgb2hsi(X);
            case 'gray'
                Y = rgb2gray(X);
            case 'hcm'
                Y = Bim_rgb2hcm(X);
            case 'pca'
                Y = Bim_rgb2pca(X);
            otherwise
                error('Bim_colorconv: conversion to %s is not supported.',s2);
        end
    case 'xyz'
        switch s2
            case 'rgb'
                Y = vl_xyz2rgb(X);
            case 'xyz'
                Y = X;
            case 'lab'
                Y = vl_xyz2lab(X);
            case 'luv'
                Y = vl_xyz2luv(X);
            case 'hsv'
                Z = vl_xyz2rgb(X);
                Y = rgb2hsv(Z);
            case 'hsi'
                Z = vl_xyz2rgb(X);
                Y = Bim_rgb2hsi(Z);
            case 'gray'
                Z = vl_xyz2rgb(X);
                Y = rgb2gray(Z);
            case 'hcm'
                Z = vl_xyz2rgb(X);
                Y = Bim_rgb2hcm(Z);
            case 'pca'
                Z = vl_xyz2rgb(X);
                Y = Bim_rgb2pca(Z);
            otherwise
                error('Bim_colorconv: conversion to %s is not supported.',s2);
        end
    case 'hsi'
        switch s2
            case 'rgb'
                Y = Bim_hsi2rgb(X);
            case 'xyz'
                Z = Bim_hsi2rgb(X);
                Y = vl_rgb2xyz(Z);
            case 'lab'
                Z = Bim_hsi2rgb(X);
                W = vl_rgb2xyz(Z);
                Y = vl_xyz2lab(W);
            case 'luv'
                Z = Bim_hsi2rgb(X);
                W = vl_rgb2xyz(Z);
                Y = vl_xyz2lub(W);
            case 'hsv'
                Z = Bim_hsi2rgb(X);
                Y = rgb2hsv(Z);
            case 'hsi'
                Y = X;
            case 'gray'
                Z = Bim_hsi2rgb(X);
                Y = rgb2gray(Z);
            case 'hcm'
                Z = Bim_hsi2rgb(X);
                Y = Bim_rgb2hcm(Z);
            case 'pca'
                Z = Bim_hsi2rgb(X);
                Y = Bim_rgb2pca(Z);
            otherwise
                error('Bim_colorconv: conversion to %s is not supported.',s2);
        end
    otherwise
        error('Bim_colorconv: conversion from %s is not supported.',s1);
end





