% F = Bsq_trifocal(P)
%
% Toolbox: Balu
%
%    Trifocal tesonsors of a sequence.
%
%    P includes the projection matrices of n views as follows:
%    Projection Pk = P(k*3-2:k*3,:), for k=1,...,n
%   
%    T are the fundamental matrices stored as follows:
%    T(p,q,r,:) are the trifocal tensors (as 27x1 vector) between views p, 
%    q and r for p=1:n-2, q=1:p+1:n-1, r=q+1:n
%
%  Example:
%
%       P1 = rand(3,4);             % proyection matrix for view 1
%       P2 = rand(3,4);             % proyection matrix for view 2
%       P3 = rand(3,4);             % proyection matrix for view 3
%       P4 = rand(3,4);             % proyection matrix for view 4
%       P5 = rand(3,4);             % proyection matrix for view 4
%       P6 = rand(3,4);             % proyection matrix for view 4
%       P = [P1;P2;P3;P4;P5;P6];    % all projection matrices
%       T = Bsq_trifocal(P);       % fundamental matrices
%       T136 = zeros(3,3,3);
%       T136(:) = T(1,3,6,:)        % Trifocal tensors between views 1-3-6
%
%  See also Bmv_fundamental, Bmv_trifocal, Bsq_fundamental.
%
% (c) D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function T = Bsq_trifocal(P)

mg = size(P,1)/3;
T = zeros(mg,mg,mg,27);
for p=1:mg-2
    p0 = 3*p-2;
    Pp = P(p0:p0+2,:);
    for q=p+1:mg-1
        q0 = 3*q-2;
        Pq = P(q0:q0+2,:);
        for r=q+1:mg
            r0 = 3*r-2;
            Pr = P(r0:r0+2,:);
            Tpqr = Bmv_trifocal(Pp,Pq,Pr);
            T(p,q,r,:) = Tpqr(:)';
        end
    end
end


