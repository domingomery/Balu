% howis(x)
%
% Toolbox: Balu
%    How is x? it displays min, max, size, class of x.
%    If x is a structure, it display the field names.
%
% D.Mery, PUC-DCC, Jul 2009-2012
% http://dmery.ing.puc.cl
%

function howis(x)
if isstruct(x)
    disp('Structure:')
    x
else
    fprintf('%s\n',['min   = ' num2str(min(x(:)))])
    fprintf('%s\n',['max   = ' num2str(max(x(:)))])
    fprintf('%s\n',['size  = ' num2str(size(x))])
    fprintf('%s\n',['class = ' class(x)])
end
