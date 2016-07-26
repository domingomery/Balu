% Toolbox: Balu
%  Display help text for Balu Matlab Toolbox
%
% D.Mery, PUC-DCC, Jun 2010
% http://dmery.ing.puc.cl
%
function Bhelp(t)
disp('Help for Balu Matlab Toolbox - (c) GRIMA, PUC-DCC')
disp(' ')
if exist('t','var')
  help BhelpImageProcessing
  help BhelpFeatureExtraction
  help BhelpFeatureTransformation
  help BhelpFeatureSelection
  help BhelpInputOutput
  help BhelpDataSelection
  help BhelpClassification
  help BhelpFeatureAnalysis
  help BhelpClustering
  help BhelpPerformanceEvaluation
  help BhelpMultipleView
  help BhelpSequenceProcessing
  help BhelpTracking
  help BhelpMiscellaneous
else
  help Bhelpmsg
  disp(' ')
  disp('>>> For all help messages type Bhelp(1)')
 
end







