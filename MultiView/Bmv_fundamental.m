% F = Bmv_fundamental(A,B,method)
%
% Toolbox: Balu
%
% Fundamental matrix from projection matrices.
%
% F = fundamental(A,B,method) returns the 3x3 fundamental matrix from
%    3x4 projection matrices A and B according to the following
%    methods:
%
%    method = 'tensor' : uses bifocal tensors with canonic matrix
%
%    method = 'tensor0': uses bifocal tensors without canonic matrix
%
%    method = 'pseudo' : uses pseudoinverse matrix of A (deault)
%
%    If no method is given, 'pseudo' will be assumed as default.
%
%    Both methods can be found in:
%
%    R. Hartley and A. Zisserman. Multiple View Geometry in Computer
%    Vision. Cambridge University Press, 2000.
%
%    Example:
%       A = rand(3,4);                      % Projection matrix for view 1
%       B = rand(3,4);                      % Projection matrix for view 2
%       M = [1 2 3 1]';                     % 3D point (X=1,Y=2,Z=3)
%       m1 = A*M;m1 = m1/m1(3);             % projection point in view 1
%       m2 = B*M;m2 = m2/m2(3);             % projection point in view 2
%       F = Bmv_fundamental(A,B,'tensor');  % Fundamental matrix using tensors
%       m2'*F*m1                            % epipolar constraint must be zero
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function F = Bmv_fundamental(A,B,method)

if ~exist('method','var')
    method = 'pseudo';
end

switch lower(method)
    
    case 'tensor'
        H = [A;rand(1,4)];
        while abs(det(H)) < 0.001
            det(H)
            H = [A;rand(1,4)];
        end
        Bs = B/H;  %Bs = B*inv(H);
        b = [Bs;Bs];
        F = zeros(3,3);
        for i=1:3
            for j=1:3
                F(i,j) = b(i+1,j)*b(i+2,4)-b(i+2,j)*b(i+1,4);
            end;
        end
        
    case 'tensor0'
        
        F = zeros(3,3);
        for i=1:3;
            sb = B;
            sb(i,:) = [];
            for j=1:3
                sa = A;
                sa(j,:) = [];
                F(i,j) = (-1)^(i+j)*det([sa; sb]);
            end;
        end
        
    otherwise
        C1 = null(A);
        Am = A'/(A*A'); % Am = A'*inv(A*A');
        F  = Bmv_antisimetric(B*C1)*B*Am;
        
        %    otherwise
        %        error('Bfundamental: method does not exist.');
end

F = F/norm(F);