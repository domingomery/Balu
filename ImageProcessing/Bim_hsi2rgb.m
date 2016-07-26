function RGB = Bim_hsi2rgb(HSI)

h = HSI(:,:,1);
s = HSI(:,:,2);
i = HSI(:,:,3);

h = 2 * pi * h;

tst1 = (0 < h) & (h <= (2 * pi/3));
tst2 = ((2 * pi/3) < h) & (h <= (4 * pi/3));
tst3 = ((4 * pi/3) < h) & (h <= (2 * pi));

f1 = (1 - s)/3;

den = cos(pi/3 - h);
den = (den == 0) * eps + (den ~= 0) .* den;

f2 = (1 + ((s .* cos(h)) ./ den )) /3;

r1 = tst1 .* f2;
g1 = tst1 .* (1 - (f1 + f2));
b1= tst1 .* f1;

den = cos(pi/3 - (h - 2 * pi/3));
den = (den == 0) * eps + (den ~= 0) .* den;

f2= (1 + ((s .* cos(h - (2 * pi/3))) ./ den)) /3;

r2 = tst2 .* f1;
g2 = tst2 .* f2;
b2 = tst2 .* (1 - (f1 + f2));

den = cos(pi/3 - (h - 4 * pi/3));
den = (den == 0) * eps + (den ~= 0) .* den;

f2= (1 + ((s .* cos(h - (4 * pi/3))) ./ den)) /3;

g3 = tst3 .* f1;
b3 = tst3 .* f2;
r3 = tst3 .* (1 - (f1 + f2));

r = r1 + r2 + r3;
g = g1 + g2 + g3;
b = b1 + b2 + b3;

r = round(255 * 3 * i .* r);
g = round(255 * 3 * i .* g);
b = round(255 * 3 * i .* b);

RGB = zeros(size(HSI));
RGB(:,:,1) = r;
RGB(:,:,2) = g;
RGB(:,:,3) = b;
