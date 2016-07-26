% F = Bmv_trifocal(A,B,C)
%
% Toolbox: Balu
%
%    Trifocal tensors.
%
%    T = trifocal(A,B,method) returns the trifocal tensors from
%    3x4 projection matrices A, B and C according of thre
%    views. T is a 3x3x3 array
%
%    The method can be found in:
%
%    R. Hartley and A. Zisserman. Multiple View Geometry in Computer
%    Vision. Cambridge University Press, 2000.
%
%    Example:
%       A = rand(3,4);           % Projection matrix for view 1
%       B = rand(3,4);           % Projection matrix for view 2
%       C = rand(3,4);           % Projection matrix for view 3
%       T = Bmv_trifocal(A,B,C)  % Trifocal tensors
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function T = Bmv_trifocal(A,B,C)

H = [A;rand(1,4)];
while abs(det(H)) < 0.001
    det(H)
    H = [A;rand(1,4)];
    if det(H) < 0.001
        disp('warning: determinant of H in trifocal estimation is zero')
    end
end
b = B/H;  %b = B*inv(H);
c = C/H;  %c = C*inv(H);
T = zeros(3,3,3);
for i=1:3;
    for j=1:3;
        for k=1:3
            T(i,j,k) = b(j,i)*c(k,4)-b(j,4)*c(k,i);
        end; 
    end; 
end
