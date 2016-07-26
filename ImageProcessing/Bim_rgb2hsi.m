function HSI = Bim_rgb2hsi(RGB)

r = double(RGB(:,:,1))/255;
g = double(RGB(:,:,2))/255;
b = double(RGB(:,:,3))/255;

i = eps + (r + g + b)/3; % Kludged to avoid problems in saturation
i = i / max(max(i));

num = .5 * ((r - g) + (r - b));
den = eps + sqrt((r - g) .* (r - g) + (r - b) .* (g - b)); %Kludged

if den == 0
   h = 0; % is this the hue we want? It doesn't matter since s = 0?
else
   h = acos(num ./ den);
end

tst = (b ./ i) > (g ./ i);
h = (2 * pi - h) .* tst + h .* (1 - tst);
h = h/(2 * pi);

s = 1 - min(r,min(g,b)) ./ i;

HSI = zeros(size(RGB));
HSI(:,:,1) = h;
HSI(:,:,2) = s;
HSI(:,:,3) = i;

