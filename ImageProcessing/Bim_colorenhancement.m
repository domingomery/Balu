function J = Bim_colorenhancement(I)

HSI = Bim_rgb2hsi(I);

HSI(:,:,3) = histeq(HSI(:,:,3));

J = Bim_hsi2rgb(HSI);

