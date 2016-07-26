% [X,Xn,options] = Bfx_lbp(I,R,options)
% [X,Xn,options] = Bfx_lbp(I,options)
% [X,Xn] = Bfx_lbp(I,R,options)
% [X,Xn] = Bfx_lbp(I,options)
%
% Toolbox: Balu
%    Local Binary Patterns features
%
%    X is the features vector, Xn is the list of feature names (see Example
%    to see how it works).
%
%    It calculates the LBP over the a regular grid of patches. The function
%    uses Heikkila & Ahonen (see http://www.cse.oulu.fi/MVG/Research/LBP).
%
%    It returns a matrix of uniform lbp82 descriptors for I, made by
%    concatenating histograms of each grid cell in the image.
%    Grid size is options.hdiv * options.vdiv
%
%    R is a binary image or empty. If R is given the lbp will be computed
%    the corresponding pixles R==0 in image I will be set to 0.
%
%     Output:
%     X is a matrix of size ((hdiv*vdiv) x 59), each row has a
%         histogram corresponding to a grid cell. We use 59 bins.
%     options.x of size hdiv*vdiv is the x coordinates of center of ith grid cell
%     options.y of size hdiv*vdiv is the y coordinates of center of ith grid cell
%     Both coordinates are calculated as if image was a square of side length 1.
%
%     References:
%     Ojala, T.; Pietikainen, M. & Maenpaa, T. Multiresolution gray-scale
%     and rotation invariant texture classification with local binary
%     patterns. IEEE Transactions on Pattern Analysis and Machine
%     Intelligence, 2002, 24, 971-987.
%
%     Mu, Y. et al (2008): Discriminative Local Binary Patterns for Human
%     Detection in Personal Album. CVPR-2008.
%
%     Example 1:
%      options.vdiv = 1;                  % one vertical divition
%      options.hdiv = 1;                  % one horizontal divition
%      options.semantic = 0;              % classic LBP
%      options.samples  = 8;              % number of neighbor samples
%      options.mappingtype = 'u2';        % uniform LBP
%      I = imread('testimg1.jpg');        % input image
%      J = I(120:219,120:239,2);          % region of interest (green)
%      figure(1);imshow(J,[])             % image to be analyzed
%      [X,Xn] = Bfx_lbp(J,[],options);    % LBP features
%      figure(2);bar(X)                   % histogram
%
%     Example 2:
%      options.vdiv = 1;                  % one vertical divition
%      options.hdiv = 1;                  % one horizontal divition
%      options.semantic = 0;              % classic LBP
%      options.samples  = 8;              % number of neighbor samples
%      options.mappingtype = 'ri';        % rotation-invariant LBP
%      I = imread('testimg1.jpg');        % input image
%      J = I(120:219,120:239,2);          % region of interest (green)
%      figure(1);imshow(J,[])             % image to be analyzed
%      [X,Xn] = Bfx_lbp(J,[],options);    % LBP features
%      bar(X)                             % histogram
%
%     Example 3:
%      options.vdiv = 1;                  % one vertical divition
%      options.hdiv = 1;                  % one horizontal divition
%      options.semantic = 1;              % semantic LBP
%      options.samples = 8;               % number of neighbor samples
%      options.sk      = 0.5;             % angle sampling
%      I = imread('testimg1.jpg');        % input image
%      J = I(120:219,120:239,2);          % region of interest (green)
%      [X,Xn] = Bfx_lbp(J,[],options);    % semantic LBP features
%      bar(X)                             % histogram
%
%     Example 4:
%      options.vdiv = 1;                  % one vertical divition
%      options.hdiv = 1;                  % one horizontal divition
%      options.semantic = 1;              % semantic LBP
%      options.samples = 16;              % number of neighbor samples
%      options.sk      = 0.5;             % angle sampling
%      I = imread('testimg1.jpg');        % input image
%      J = I(120:219,120:239,2);          % region of interest (green)
%      [X,Xn] = Bfx_lbp(J,[],options);    % semantic LBP features
%      bar(X)                             % histogram
%
%     Example 5:
%      options.vdiv = 1;                  % one vertical divition
%      options.hdiv = 1;                  % one horizontal divition
%      options.semantic = 1;              % semantic LBP
%      options.samples = 8;               % number of neighbor samples
%      options.sk      = 0.25;            % angle sampling
%      options.weight  = 9;               % angle sampling
%      I = imread('testimg1.jpg');        % input image
%      J = I(120:219,120:239,2);          % region of interest (green)
%      [X,Xn] = Bfx_lbp(J,[],options);    % weighted LBP features
%      bar(X)                             % histogram
%   See also Bfx_gabor, Bfx_clp, Bfx_fourier, Bfx_dct.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl
%
% James Kapaldo updated this version for Matlab2014b
%
function [X,Xn,options] = Bfx_lbp(I,R,options)

if nargin==2;
    options = R;
    R = ones(size(I));
end

vdiv = options.vdiv;
hdiv = options.hdiv;

if ~isfield(options,'show')
    options.show = 0;
end

if ~isfield(options,'normalize')
    options.normalize = 0;
end



if options.show == 1
    disp('--- extracting local binary patterns features...');
end

if ~isfield(options,'samples')
    options.samples = 8;
end

if ~isfield(options,'integral')
    options.integral = 0;
end



if ~isfield(options,'radius')
    options.radius = log(options.samples)/log(2)-1;
end

if ~isfield(options,'semantic')
    options.semantic = 0;
end

if ~isfield(options,'weight')
    options.weight = 0;
end

LBPst = 'LBP';

if options.semantic>0
    if ~isfield(options,'sk')
        options.sk = 1;
    end
    mapping = getsmapping(options.samples,options.sk);
    LBPst = ['s' LBPst];
    st='8x8';
else
    %    mapping = getmapping(8,'u2');
    if ~isfield(options,'mappingtype')
        options.mappingtype = 'u2';
    end
    st = sprintf('%d,%s',options.samples,options.mappingtype);
    mapping = getmapping(options.samples,options.mappingtype);
end
% get lbp image




if ~isempty(R);
    I(R==0) = 0;
end
code_img = lbp(I,options.radius,options.samples,mapping,'');
[n1,n2]  = size(code_img);
[N,M]    = size(I);
Ilbp     = zeros(size(I));
i1 = round((N-n1)/2);
j1 = round((M-n2)/2);
Ilbp(i1+1:i1+n1,j1+1:j1+n2) = code_img;
options.Ilbp = Ilbp;

if options.integral == 1
    options.Hx = Bim_inthist(Ilbp+1,options.maxD);
end



ylen     = round(n1/vdiv);
xlen     = round(n2/hdiv);
% split image into blocks (saved as columns)
grid_img = im2col(code_img,[ylen, xlen], 'distinct');
if options.weight>0
    LBPst = ['w' LBPst];
    mt  = 2*options.radius-1;
    mt2 = mt^2;
    Id = double(I);
    switch options.weight
        case 1
            W = abs(conv2(Id,ones(mt,mt)/mt2,'same')-Id);
        case 2
            W = (abs(conv2(Id,ones(mt,mt)/mt2,'same')-Id))./(Id+1);
        case 3
            W = abs(medfilt2(Id,[mt mt])-Id);
        case 4
            W = abs(medfilt2(Id,[mt mt])-Id)./(Id+1);
        case 5
            W = abs(ordfilt2(Id,mt2,ones(mt,mt))-Id);
        case 6
            W = abs(ordfilt2(Id,mt2,ones(mt,mt))-Id)./(Id+1);
        case 7
            Id = conv2(Id,ones(mt,mt)/mt2,'same');
            W = abs(ordfilt2(Id,mt2,ones(mt,mt))-Id)./(Id+1);
        case 8
            Id = medfilt2(Id,[mt mt]);
            W = abs(ordfilt2(Id,mt2,ones(mt,mt))-Id)./(Id+1);
        case 9
            Id = medfilt2(Id,[mt mt]);
            W = abs(ordfilt2(Id,mt2-1,ones(mt,mt))-Id)./(Id+1);
        otherwise
            error('Bfx_lbp does not recognice options.weight = %d.',options.weight);
    end
    W = W(mt+1:end-mt,mt+1:end-mt);
    grid_W = im2col(W,[ylen, xlen], 'distinct');
    nwi = mapping.num;
    nwj = size(grid_W,2);
    nwk = size(grid_W,1);
    desc = zeros(nwi,nwj);
    for j=1:nwj
        x = grid_img(:,j)+1;
        y = grid_W(:,j);
        d = zeros(nwi,1);
        for k=1:nwk
            d(x(k))=d(x(k))+y(k);
            % d(x(k))=d(x(k))+1; % normal LBP each LBP has equal weight
        end
        desc(:,j) = d;
    end
    
else
    desc = hist(double(grid_img), 0:mapping.num-1);
    % calculate coordinates of descriptors as if I was square w/ side=1
end
dx = 1.0/hdiv;
dy = 1.0/vdiv;
x = dx/2.0: dx :1.0-dx/2.0;
y = dy/2.0: dy :1.0-dy/2.0;
options.x = x;
options.y = y;
if hdiv*vdiv>1
    D = desc';
else
    D = desc;
end
[M,N] = size(D);

Xn = char(zeros(N*M,24));
X  = zeros(1,N*M);
k=0;
for i=1:M
    for j=1:N
        k = k+1;
        s = sprintf('%s(%d,%d)[%s]                       ',LBPst,i,j,st);
        Xn(k,:) = s(1:24);
        X(k) = D(i,j);
    end
end

if options.normalize
    X = X/sum(X);
end


end


%GETMAPPING returns a structure containing a mapping table for LBP codes.
%  MAPPING = GETMAPPING(SAMPLES,MAPPINGTYPE) returns a
%  structure containing a mapping table for
%  LBP codes in a neighbourhood of SAMPLES sampling
%  points. Possible values for MAPPINGTYPE are
%       'u2'   for uniform LBP
%       'ri'   for rotation-invariant LBP
%       'riu2' for uniform rotation-invariant LBP.
%
%  Example:
%       I=imread('rice.tif');
%       MAPPING=getmapping(16,'riu2');
%       LBPHIST=lbp(I,2,16,MAPPING,'hist');
%  Now LBPHIST contains a rotation-invariant uniform LBP
%  histogram in a (16,2) neighbourhood.
%

function mapping = getmapping(samples,mappingtype)
% Version 0.1.1
% Authors: Marko Heikkila and Timo Ahonen

% Changelog
% 0.1.1 Changed output to be a structure
% Fixed a bug causing out of memory errors when generating rotation
% invariant mappings with high number of sampling points.
% Lauge Sorensen is acknowledged for spotting this problem.



table = 0:2^samples-1;
newMax  = 0; %number of patterns in the resulting LBP code
index   = 0;

%vr2014b = or(strcmp(version('-release'),'2014b'),strcmp(version('-release'),'2014a'));
%if vr2014b
    switch samples
        case 8
            sampleType = 'uint8';
        case 16
            sampleType = 'uint16';
        otherwise
    end
%else
%    sampleType = samples;
%end

if strcmp(mappingtype,'u2') %Uniform 2
    newMax = samples*(samples-1) + 3;
    for i = 0:2^samples-1
        j = bitset(bitshift(i,1,sampleType),1,bitget(i,samples)); %rotate left
        numt = sum(bitget(bitxor(i,j),1:samples)); %number of 1->0 and
        %0->1 transitions
        %in binary string
        %x is equal to the
        %number of 1-bits in
        %XOR(x,Rotate left(x))
        if numt <= 2
            table(i+1) = index;
            index = index + 1;
        else
            table(i+1) = newMax - 1;
        end
    end
end

if strcmp(mappingtype,'ri') %Rotation invariant
    tmpMap = zeros(2^samples,1) - 1;
    for i = 0:2^samples-1
        rm = i;
        r  = i;
        for j = 1:samples-1
            r = bitset(bitshift(r,1,sampleType),1,bitget(r,samples)); %rotate
            %left
            if r < rm
                rm = r;
            end
        end
        if tmpMap(rm+1) < 0
            tmpMap(rm+1) = newMax;
            newMax = newMax + 1;
        end
        table(i+1) = tmpMap(rm+1);
    end
end

if strcmp(mappingtype,'riu2') %Uniform & Rotation invariant
    newMax = samples + 2;
    for i = 0:2^samples - 1
        j = bitset(bitshift(i,1,sampleType),1,bitget(i,samples)); %rotate left
        numt = sum(bitget(bitxor(i,j),1:samples));
        if numt <= 2
            table(i+1) = sum(bitget(i,1:samples));
        else
            table(i+1) = samples+1;
        end
    end
end

mapping.table=table;
mapping.samples=samples;
mapping.num=newMax;
end


%  LBP returns the local binary pattern image or LBP histogram of an image.
%  J = LBP(I,R,N,MAPPING,MODE) returns either a local binary pattern
%  coded image or the local binary pattern histogram of an intensity
%  image I. The LBP codes are computed using N sampling points on a
%  circle of radius R and using mapping table defined by MAPPING.
%  See the getmapping function for different mappings and use 0 for
%  no mapping. Possible values for MODE are
%       'h' or 'hist'  to get a histogram of LBP codes
%       'nh'           to get a normalized histogram
%  Otherwise an LBP code image is returned.
%
%  J = LBP(I) returns the original (basic) LBP histogram of image I
%
%  J = LBP(I,SP,MAPPING,MODE) computes the LBP codes using n sampling
%  points defined in (n * 2) matrix SP. The sampling points should be
%  defined around the origin (coordinates (0,0)).
%
%  Examples
%  --------
%       I=imread('rice.png');
%       mapping=getmapping(8,'u2');
%       H1=LBP(I,1,8,mapping,'h'); %LBP histogram in (8,1) neighborhood
%                                  %using uniform patterns
%       subplot(2,1,1),stem(H1);
%
%       H2=LBP(I);
%       subplot(2,1,2),stem(H2);
%
%       SP=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
%       I2=LBP(I,SP,0,'i'); %LBP code image using sampling points in SP
%                           %and no mapping. Now H2 is equal to histogram
%                           %of I2.

function result = lbp(varargin) % image,radius,neighbors,mapping,mode)
% Version 0.3.2
% Authors: Marko Heikkila and Timo Ahonen

% Changelog
% Version 0.3.2: A bug fix to enable using mappings together with a
% predefined spoints array
% Version 0.3.1: Changed MAPPING input to be a struct containing the mapping
% table and the number of bins to make the function run faster with high number
% of sampling points. Lauge Sorensen is acknowledged for spotting this problem.

% Check number of input arguments.
error(nargchk(1,5,nargin));

image=varargin{1};
d_image=double(image);

if nargin==1
    spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
    neighbors=8;
    mapping=0;
    mode='h';
end

if (nargin == 2) && (length(varargin{2}) == 1)
    error('Input arguments');
end

if (nargin > 2) && (length(varargin{2}) == 1)
    radius=varargin{2};
    neighbors=varargin{3};
    
    spoints=zeros(neighbors,2);
    
    % Angle step.
    a = 2*pi/neighbors;
    
    for i = 1:neighbors
        spoints(i,1) = -radius*sin((i-1)*a);
        spoints(i,2) = radius*cos((i-1)*a);
    end
    
    if(nargin >= 4)
        mapping=varargin{4};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 5)
        mode=varargin{5};
    else
        mode='h';
    end
end

if (nargin > 1) && (length(varargin{2}) > 1)
    spoints=varargin{2};
    neighbors=size(spoints,1);
    
    if(nargin >= 3)
        mapping=varargin{3};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 4)
        mode=varargin{4};
    else
        mode='h';
    end
end

% Determine the dimensions of the input image.
[ysize xsize] = size(image);



miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;

% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0));
origx=1-floor(min(minx,0));

% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(xsize < bsizex || ysize < bsizey)
    error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;

% Fill the center pixel matrix C.
C = image(origy:origy+dy,origx:origx+dx);
d_C = double(C);

bins = 2^neighbors;

% Initialize the result matrix with zeros.
result=zeros(dy+1,dx+1);

%Compute the LBP code image

for i = 1:neighbors
    y = spoints(i,1)+origy;
    x = spoints(i,2)+origx;
    % Calculate floors, ceils and rounds for the x and y.
    fy = floor(y); cy = ceil(y); ry = round(y);
    fx = floor(x); cx = ceil(x); rx = round(x);
    % Check if interpolation is needed.
    if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
        % Interpolation is not needed, use original datatypes
        N = image(ry:ry+dy,rx:rx+dx);
        D = N >= C;
    else
        % Interpolation needed, use double type images
        ty = y - fy;
        tx = x - fx;
        
        % Calculate the interpolation weights.
        w1 = (1 - tx) * (1 - ty);
        w2 =      tx  * (1 - ty);
        w3 = (1 - tx) *      ty ;
        w4 =      tx  *      ty ;
        % Compute interpolated pixel values
        N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
            w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
        D = N >= d_C;
    end
    % Update the result matrix.
    v = 2^(i-1);
    result = result + v*D;
end

%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    for i = 1:size(result,1)
        for j = 1:size(result,2)
            result(i,j) = mapping.table(result(i,j)+1);
        end
    end
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    result=hist(result(:),0:(bins-1));
    if (strcmp(mode,'nh'))
        result=result/sum(result);
    end
else
    %Otherwise return a matrix of unsigned integers
    if ((bins-1)<=intmax('uint8'))
        result=uint8(result);
    elseif ((bins-1)<=intmax('uint16'))
        result=uint16(result);
    else
        result=uint32(result);
    end
end

end



% Mapping for sLBB
function mapping = getsmapping(N,sk)


vr2014b = or(strcmp(version('-release'),'2014b'),strcmp(version('-release'),'2014a'));
if vr2014b
    switch N
        case 8
            sampleType = 'uint8';
        case 16
            sampleType = 'uint16';
        otherwise
    end
else
    sampleType = samples;
end


M = 2^N;
samples = N;
len = zeros(M,1);
ang = zeros(M,1);
for x=0:(M-1)
    k = x+1;
    j = bitset(bitshift(x,1,sampleType),1,bitget(x,samples));
    numt = sum(bitget(bitxor(x,j),1:samples));
    c = numt;
    if c>2
        len(k)=-1;
        ang(k)=-1;
    else
        s = bitget(x,1:samples);
        len(k) = sum(s);
        if c==0
            ang(k)=0;
        else
            r = 0;
            while (s(1)~=0) || (s(N)~=1)
                s = [s(2:N) s(1)];
                r = r+1;
            end
            ii = find(s==1);
            a = mean(ii)+r;
            if a>N
                a=a-N;
            end
            ang(k) = round(sk*a)-1;
        end
    end
    % fprintf('%4d: %s (%d,%d)\n',x,dec2bin(x,N),len(k),ang(k)); pause
end

Ma = max(ang)+1;
map = len*Ma+ang-Ma+1;
n = max(map)+1;
map(ang==-1) = n;
map(1) = 0;



mapping.table   = map';
mapping.samples = N;
mapping.num     = n+1;


end




