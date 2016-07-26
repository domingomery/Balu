% D.Mery
% Ago-2009
% Segmentacion de fallas hipoteticas
% input : X imagen de rayos X de entrada
% output: Y imagen binaria de salida con las fallas hipoteticas
% output: feat matriz de features extraidas ver ultimas lineas del codigo

function [Y,feat] = Bim_segdefects(X)

%
%Inicializacion de parametros:
%

%ajustes de contraste para imadjust
pb   = [30 0190 1.5];% 1: gris minimo
% 2: gris maximo
% 3: gamma

%deteccion de bordes
pe  = [1 11 0 15 10]; % 1: 1=>LoG, 2=>Canny, 3=>LoGX
% 2: Tamanho de la mascara
% 3: Umbral de la deteccion de bordes
% 4: Mascara del gradiente
% 5: Umbral del gradiente

%filtro mediana
pm = [15 3];          % 1: Tamanho de la mascara
% 2: umbral
%pre filtro
%pf  = [1  3];         %1: 0=>no hay, 1=>Gauss, 2=>Median, 3=>promedio
%2: Mascara

%clasificacion
% pk = [-0.0056 0.45];

%Otros:
Amin = 7;             % area minima en pixeles
Amax = 250;           % area maxima en pixeles
Bmax = 230;           % maximo valor de intensidad permitido (si es muy grande es un agujero)
Bmin = 30;            % minimo valor de intensidad permitido
Dmax = -0.7;          % maximo valor de la segunda derivada (valores negativos-> fallas claras)
Ksth = 2;             % umbral para el contraste Ks
Fth  = 0.65;          % umbral para el CLP
RAth = 0.05;          % umbral para la razon de areas


%
% ajuste de contraste de la imagen
%

xmin = pb(1);
xmax = pb(2);
B = abs(double(X)-xmin)/(xmax-xmin)*255; % expansion lineal
B = imadjust(B/256,[],[],pb(3))*255+1;   % ajuste gamma

%
% deteccion de fallas con mediana
%
Y = (B - medfilt2(B,[pm(1) pm(1)]))>pm(2);

%
% preprocesamiento y definiciones
%
%F         = conv2(B,fspecial('gaussian',pf(2),pf(2)/8.5),'same');  % pasa bajos
Ea        = edge(B,'log',pe(3),pe(2)/8.5);                         % deteccion de bordes
G         = Bim_d1(B,pe(4));          % gradiente de la imagen
E         = not(Ea|G>pe(5));          % ensanche de bordes con gradiente alto
G         = bwareaopen(E,Amin+1,4);   % eliminacion de regiones pequenhas
Go        = imclearborder(G,4);       % eliminacion de regiones que tocan los limites de la imagen
L         = bwlabel(Go,4);            % deteccion de regiones cerradas
n         = max(L(:));                % numero de regiones
X2        = Bim_d2(B);                % segunda derivada
[N,M]     = size(B);                  % dimensiones de la imagen
Dete      = zeros(N,M);               % imagen donde se almacenaran las fallas hipoteticas detectadas
feat      = zeros(10000,11);
it        = 0;


stats = regionprops(L,'PixelIdxList','PixelList','Centroid');

%
% extraccion de caracteristicas de cada region cerrada encontrada en la imagen L
% y segmentacion
%
for i = 1:n
    j      = stats(i).PixelIdxList;   % busca los pixeles con label i
    D      = mean(X2(j));             % valor medio de la segunda derivada de la region
    if ((D<Dmax))                     % es la segunda aderivada ok?
        G      = mean(B(j));          % valor medio de grises de la region
        if ((G<Bmax) && (G>Bmin))     % es el tono de gris
            A      = length(j);       % area
            RA     = sum(Y(j))/A;     % pixeles detectados por el filtro mediana / area (Mery, 2006: High contrast pixels)
            if (A<Amax)&&(RA>RAth)    % area minima y razon de areas
                
                REG    = zeros(N,M);                           % se inicializa una imagen REG vacia
                REG(j) = 1;                                    % se marca en REG los pixeles de la region i
                i_max       = max(stats(i).PixelList(:,2))+1;  % maximo en la direccion i
                i_min       = min(stats(i).PixelList(:,2))-1;  % minimo en la direccion i
                j_max       = max(stats(i).PixelList(:,1))+1;  % maximo de en la direccion j
                j_min      =  min(stats(i).PixelList(:,1))-1;  % minimo de en la direccion j
                
                h           = i_max-i_min+1;                   % alto
                b           = j_max-j_min+1;                   % ancho
                i_m         = stats(i).Centroid(2);            % centro de gravedad i
                j_m         = stats(i).Centroid(1);            % centro de gravedad j
                
                % Definicion de una ventana que incluya la region y su entorno
                imf = fix(i_m+0.5);
                jmf = fix(j_m+0.5);
                Ia = imf-h; Ie = imf+h;
                Ja = jmf-b; Je = jmf+b;
                if ((Ia>1)&&(Ie<N)&&(Ja>1)&&(Je<M)) % la ventana no debe salir de la imagen
                    
                    Gb    = B(Ia:Ie,Ja:Je);                    % ventana
                    op.show = 0;
                    op.neighbor = 1;
                    op.param = 1.5;
                    X = Bfx_contrast(Gb,ones(size(Gb)),op);
                    Ks = X(4);
                    op.ng = 32;
                    X = Bfx_clp_old(Gb,ones(size(Gb)),op);
                    
                    FTKn2 = X(6);
                    
                    
                    if (FTKn2>Fth)&&(Ks>Ksth)
                        Dete = or(Dete,REG);
                        it = it+1;
                        feat(it,:) = [i_m j_m sqrt(A) h b A G D FTKn2 Ks RA];
                    end
                end
            end
        end
    end
end
Y = Dete;
feat = feat(1:it,:);
