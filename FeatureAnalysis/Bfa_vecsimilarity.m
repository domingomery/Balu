% [rk,j] = Bvecsimilarity(vq,v)
%
% Toolbox: Balu
%
%     Normalized scalar product (cosine angle).
%
%     vq: query vector
%     v : family of vectors
%     This function serach the minimal distance between vq and all vector
%     in v.
%     rk are the sorted scalar product.
%     j are the sorted indices.
%     
%  Example:
%     See examples from Bsq_vgoogle.
%
%  Reference:
%     Sivic & Zisserman: Efficient visual search for videos cast as text
%     retrieval. 31(4):591-606. PAMI 2009.
%
%  See also Bseqvgoogle.
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl


function [rk,j] = Bfa_vecsimilarity(vq,v)

N = size(v,1);
nvq = norm(vq);
if nvq>0
    s = zeros(N,1);
    for d=1:N
        vd = v(d,:)';
        nvd = norm(vd);
        if nvd>0
            s(d) = abs(vq'*vd)/nvq/nvd;
        end
    end
    [rk,j] = sort(s,'descend');
else
    error('Bfa_vecsimilarity: Norm of vq is zero...');
end
