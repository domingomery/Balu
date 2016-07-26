%function st = num2fixstr(i,d)
%
% This function converts an integer i in a string
% st with a fixed number of characters (with '0'
% filled from left.
%
% Example:
%     st = num2fixstr(3,4) % returns '0003'
%
% D.Mery, Aug-2012
% http://dmery.ing.puc.cl
%

function st = num2fixstr(i,d)
s = [char('0'*ones(1,d)) num2str(i)];
st = s(end-d+1:end);