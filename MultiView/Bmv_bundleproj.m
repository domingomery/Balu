% [Xs,Ps,xs] = Bmv_bundleproj(x)
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
%       Ps estimation of projection matrices as 3m x 4 matrix
%       xs reprojected points as m x n x 2 matrix
%
%       lambda * [x y 1]' = P*[X Y Z 1]'  > Projective projection
%
% Reference: 
%    Hartley & Zisserman: Multiple View Geometry, 2nd Edition, p. 444
%
%    Example:
%       n = 20; m = 6;                  % 20 points in 6 views
%       P1 = rand(3,4);                 % Projection matrix for view 1
%       P2 = rand(3,4);                 % Projection matrix for view 2
%       P3 = rand(3,4);                 % Projection matrix for view 3
%       P4 = rand(3,4);                 % Projection matrix for view 4 
%       P5 = rand(3,4);                 % Projection matrix for view 5
%       P6 = rand(3,4);                 % Projection matrix for view 6
%       M = [rand(3,n);ones(1,n)];      % n 3D points
%       m1 = P1*M;                      % projection points in view 1
%       m2 = P2*M;                      % projection points in view 2
%       m3 = P3*M;                      % projection points in view 3
%       m4 = P4*M;                      % projection points in view 4
%       m5 = P5*M;                      % projection points in view 5
%       m6 = P6*M;                      % projection points in view 6
%       x = zeros(3,n,m);
%       x(:,:,1) = m1./(ones(3,1)*m1(3,:));
%       x(:,:,2) = m2./(ones(3,1)*m2(3,:));
%       x(:,:,3) = m3./(ones(3,1)*m3(3,:));
%       x(:,:,4) = m4./(ones(3,1)*m4(3,:));
%       x(:,:,5) = m5./(ones(3,1)*m5(3,:));
%       x(:,:,6) = m6./(ones(3,1)*m6(3,:));
%       [Xs,Ps,xs] = Bmv_bundleproj(x); % Bundle Adjustment
%       e = xs-x                        % reprojection error
%
% See also Bmv_bundleafin.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [Xs,Ps,xs] = Bmv_bundleproj(x)

% Normalization

xn = x;
xc = x(1:2,:,:);
mx = max(xc(:));
xn(1:2,:,:) = x(1:2,:,:)/mx;


K = [mx 0 0; 0 mx 0; 0 0 1];


% Definitions
[p,n,m] = size(x);

lambda = ones(m,n);

s0 = Inf;

ok = 0;
t = 0;
while not(ok)
    t = t+1;
    
    for i=1:m
        lambda(i,:) = lambda(i,:)/max(lambda(i,:));
    end
    
    for j=1:n
        lambda(:,j) = lambda(:,j)/max(lambda(:,j));
    end
    
    % Construction of measurement matrix W
    W = zeros(3*m,n);
    for i=1:m
        for j=1:n
            W(i*3-2:i*3,j) = lambda(i,j)*xn(:,j,i);
        end
    end
    [U,S,V] = svd(W);

    % Estimation of Ps and Xs
    Uf = U(:,1:4);
    Vf = V(:,1:4);
    Df = S(1:4,1:4);

    Ps = Uf*Df;
    Xs = Vf';

    xs = zeros(3,n,m); % re-projected 2D points
    for i=1:m
        for j=1:n
            xs(:,j,i)   = Ps(i*3-2:3*i,:)*Xs(:,j);
            lambda(i,j) = xs(3,j,i);
            xs(:,j,i)   = K*xs(:,j,i)/xs(3,j,i);
        end
    end

    % Geometric error
    s = 0;
    d = zeros(2,1);
    for i=1:m
        for j=1:n
            d(:) = abs(x(1:2,j,i)-xs(1:2,j,i));
            s = s + norm(d);
        end
    end

    s = s/n/m;
    ok = or(abs(s-s0)<1e-4,t>15);
    s0 = s;
end


for i=1:m
    Ps(i*3-2:3*i,:) = K*Ps(i*3-2:3*i,:);
end

