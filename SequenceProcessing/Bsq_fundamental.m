% F = Bsq_fundamental(P)
%
% Toolbox: Balu
%
%    Fundamental matrices of a sequence.
%
%    P includes the projection matrices of n views as follows:
%    Projection Pk = P(k*3-2:k*3,:), for k=1,...,n
%   
%    F are the fundamental matrices stored as follows:
%    F(p,q,:) is the Fundamental matrix (as 9x1 vector) between view p and
%    view q, for p=1:n-1 and q=1:p+1:n
%
%  Example:
%
%       M = [1 2 3 1]';             % 3D point (X=1,Y=2,Z=3)
%       P1 = rand(3,4);             % proyection matrix for view 1
%       P2 = rand(3,4);             % proyection matrix for view 2
%       P3 = rand(3,4);             % proyection matrix for view 3
%       P4 = rand(3,4);             % proyection matrix for view 4
%       m1 = P1*M; m1=m1/m1(3);     % proyection point in view 1
%       m2 = P2*M; m2=m2/m2(3);     % proyection point in view 2
%       m3 = P3*M; m3=m3/m3(3);     % proyection point in view 3
%       m4 = P4*M; m4=m4/m4(3);     % proyection point in view 4
%       P = [P1;P2;P3;P4];          % all projection matrices
%       F = Bsq_fundamental(P);    % fundamental matrices
%       F13 = zeros(3,3);
%       F13(:) = F(1,3,:);          % Fundamental matrix between views 1-3
%       m3'*F13*m1                  % epipolar constraint must be zero
%       F24 = zeros(3,3);
%       F24(:) = F(2,4,:);          % Fundamental matrix between views 2-4
%       m4'*F24*m2                  % epipolar constraint must be zero
%
%  See also Bmv_fundamental, Bmv_trifocal, Bsq_trifocal.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function F = Bsq_fundamental(P)

n = size(P,1)/3;

F = zeros(n,n,9);

for p=1:n-1
    i0 = 3*p-2;
    Pi = P(i0:i0+2,:);
    for q=p+1:n
        j0 = 3*q-2;
        Pj = P(j0:j0+2,:);
        Fij = Bmv_fundamental(Pi,Pj);
        F(p,q,:) = Fij(:)';
    end
end

