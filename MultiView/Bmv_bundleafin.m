% [Xs,Ps,xs] = Bmv_bundleafin(x)
%
% Toolbox: Balu
%
%    Bundle Adjustment Projective reconstruction 
%    using the factorization algorithm
%
%    input:
%       x projected 2D points as 3 x n x m matrix (homogeneous)
%         with n number of 3D points
%              m number of views
%
%    output:
%       Xs estimation of 3D points as 3 x n matrix
%       Ps estimation of affine projection matrix Ps = [Ms ts; 0 0 0 1]
%          Ms estimation of affine matrices as 2m x 3 matrix
%          ts estimation of affine vectors as 2 x m matrix
%       xs reprojected points as m x n x 2 matrix
%
%      [x y]' = M*[X Y Z]' + t > Affine projection
%
% Reference: 
%    Hartley & Zisserman: Multiple View Geometry, 2nd Edition, p. 437
%
%    Example:
%       n = 20; m = 6;                  % 20 points in 6 views
%       P1 = [rand(2,4);[0 0 0 1]];     % Projection matrix for view 1
%       P2 = [rand(2,4);[0 0 0 1]];     % Projection matrix for view 2
%       P3 = [rand(2,4);[0 0 0 1]];     % Projection matrix for view 3
%       P4 = [rand(2,4);[0 0 0 1]];     % Projection matrix for view 4
%       P5 = [rand(2,4);[0 0 0 1]];     % Projection matrix for view 5
%       P6 = [rand(2,4);[0 0 0 1]];     % Projection matrix for view 6
%       M = [rand(3,n);ones(1,n)];      % n 3D points
%       x(:,:,1) = P1*M;                % projection points in view 1
%       x(:,:,2) = P2*M;                % projection points in view 2
%       x(:,:,3) = P3*M;                % projection points in view 3
%       x(:,:,4) = P4*M;                % projection points in view 4
%       x(:,:,5) = P5*M;                % projection points in view 5
%       x(:,:,6) = P6*M;                % projection points in view 6
%       [Xs,Ps,xs] = Bmv_bundleafin(x); % Bundle Adjustment
%       e = xs-x                        % reprojection error
% 
% See also Bmv_bundleproj.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [Xs,Ps,xs] = Bmv_bundleafin(x)

nx = size(x,1);

% Definitions
[p,n,m] = size(x);
ts      = zeros(2,m);

% (i) Computation of translations:
for i=1:m
    mtsi = mean(x(:,:,i),2);
    ts(:,i) = mtsi(1:2);
end

% (ii) Centre the data:
xc = x;
c = zeros(2,1);
for i=1:m
    for j = 1:n
        c(:)      = x(1:2,j,i);
        xc(1:2,j,i) = c-ts(:,i);
    end
end

% (iii) Constructionof measurement matrix W
W = zeros(2*m,n);
for i=1:m
    for j=1:n
        W(i*2-1:i*2,j) = xc(1:2,j,i);
    end
end
[U,S,V] = svd(W);

% (iv) Estimation of Ms and Xs
Uf = U(:,1:3);
Vf = V(:,1:3);
Df = S(1:3,1:3);

Ms = Uf*Df;
Xs = Vf';

xs = ones(nx,n,m); % re-projected 2D points
for i=1:m
    for j=1:n
        xs(1:2,j,i) = Ms(i*2-1:i*2,:)*Xs(:,j)+ts(:,i);
    end
end

Ps = zeros(3*m,4);
for i=1:m
    Ps(i*3-2:i*3,:) = [Ms(i*2-1:i*2,:) ts(:,i); 0 0 0 1];
end


