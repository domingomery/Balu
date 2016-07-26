% Obtain the index number of the feature wich names contain a given string.

function [ix,fnix] = Bio_findex(fn,str,inc)

if ~exist('inc','var')
    inc = 1;
end


[N,M] = size(fn);
n = length(str);

T=ones(N,1)*str;
ix = [];
for i=1:M-n+1
   D  = sum(abs(T-fn(:,i:i+n-1))')';
   ii = find(D==0);
   if ~isempty(ii)
       ix = [ix;ii];
   end
end
if not(inc)
    ii     = (1:N)';
    ii(ix) = [];
    ix     = ii;
end
fnix = fn(ix,:);
