% sc = Bfa_sqcorrcoef(x,y)
%
% Toolbox: Balu
%    Squared-correlation coefficient between two random vector x and y
%
%    Wei, H.-L. & Billings, S. Feature Subset Selection and Ranking for
%    Data Dimensionality Reduction Pattern Analysis and Machine 
%    Intelligence, IEEE Transactions on, 2007, 29, 162-166
%
% D.Mery, PUC-DCC, Jul. 2009
% http://dmery.ing.puc.cl
%
function sc = Bfa_sqcorrcoef(x,y)

sc = ((x'*y)/norm(x)/norm(y))^2;

end