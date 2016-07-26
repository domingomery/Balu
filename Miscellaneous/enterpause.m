% Toolbox: Balu
%    enterpause display "press <Enter> to continue..." and wait for <Enter>
%    enterpause(t) means pause(t)
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl
function enterpause(t)
if ~exist('t','var')
    disp('press <Enter> to continue...')
    pause
else
    pause(t)
end