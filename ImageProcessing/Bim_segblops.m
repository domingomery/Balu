% D.Mery
% Octo-2010
% Segmentacion de blops
% input : X imagen de entrada
% output: Y imagen binaria de salida con los blops detectados
% output: feat matriz de features extraidas ver ultimas lineas del codigo

function [Y,feat] = Bim_segblops(X,pa,show)


if ~exist('show','var')
    show = 0;
end


n    = pa(1);  % mascara mediana
t    = pa(2);  % umbral para la resta original-mediana
Bmin = pa(3);  % minimo valor de gris permitido
Bmax = pa(4);  % maximo valor de gris permitido
Amin = pa(5);  % area minima
Amax = pa(6);  % area maxima
Emin = pa(7);  % solidez minima
Omin = pa(8);  % contraste minimo

S = medfilt2(X,[n n]);
J = abs(S-X)>t;  % umbralizacion del filtro mediana
K = and(J,(X<Bmax)&(X>Bmin));    % seleccion del rango de los tonos de gris
P = bwareaopen(K,Amin,4);        % eliminacion de regiones pequenhas
Q = imclearborder(P,4);          % eliminacion de regiones que tocan los limites de la imagen
L = bwlabel(Q,4);                % deteccion de regiones cerradas
n = max(L(:));                   % numero de regiones


stats = regionprops(L,'PixelIdxList','PixelList','Centroid','Solidity');

[N,M]     = size(X);             % dimensiones de la imagen
Y         = zeros(N,M);          % imagen donde se almacenaran las blops detectados
feat      = zeros(10000,4);
it        = 0;

%
% segmentacion: evaluacion para cada region i = 1...n
%

for i = 1:n
    j      = stats(i).PixelIdxList;             % busca los pixeles con label i
    e      = stats(i).Solidity;                 % solidez de la region
    G      = mean(X(j));                        % valor medio de grises
    if ((G<Bmax) && (G>Bmin) && (e>Emin))
        A      = length(j);                     % area
        if (A<Amax)
            i_m        = stats(i).Centroid(2);  % centro de gravedad i
            j_m        = stats(i).Centroid(1);  % centro de gravedad j
            i_max       = max(stats(i).PixelList(:,2))+1;  % maximo en la direccion i
            i_min       = min(stats(i).PixelList(:,2))-1;  % minimo en la direccion i
            j_max       = max(stats(i).PixelList(:,1))+1;  % maximo de en la direccion j
            j_min      =  min(stats(i).PixelList(:,1))-1;  % minimo de en la direccion j
            h           = i_max-i_min+1;                   % alto
            b           = j_max-j_min+1;                   % ancho

            % Definicion de una ventana que incluya la region y su entorno
            imf = fix(i_m+0.5);
            jmf = fix(j_m+0.5);
            Ia = imf-h; Ie = imf+h;
            Ja = jmf-b; Je = jmf+b;

            if ((Ia>1)&&(Ie<N)&&(Ja>1)&&(Je<M)) % la ventana no debe salir de la imagen
                REG    = zeros(N,M);            % se inicializa una imagen REG vacia
                REG(j) = 1;                     % se marca en REG los pixeles de la region i

                Gb    = X(Ia:Ie,Ja:Je);         % ventana
                REGb  = REG(Ia:Ie,Ja:Je);       % ventana de la region (imagen binaria)
                
                % Contraste
                REGe  = and(not(REGb),(imdilate(REGb,ones(5,5))));
                o = abs(mean(Gb(REGb==1))-mean(Gb(REGe==1)))/mean(Gb(REGb==1));

                
                
                if o>Omin
                    Y(j)       = 1;
                    it         = it+1;
                    feat(it,:) = [i_m j_m sqrt(A) o];
                end
                %figure(1)
                %imshow(Gb,[]);
                %figure(2)
                %imshow(REGe,[])
                %figure(3)
                %imshow(REGb,[])
                %enterpause
                
            end
        end
    end
end
feat = feat(1:it,:);


if show
    figure(1);imshow(X,[]);title('X: original');
    figure(2);imshow(S,[]);title('S: median');
    figure(3);imshow(J,[]);title('J = |S-X|>t');
    figure(4);imshow(P,[]);title('J without small regions');
    figure(5);Bio_edgeview(X,Y)
end