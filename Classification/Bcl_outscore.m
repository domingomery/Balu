function ds = Bcl_outscore(ds,score,options)
outp = 0;
if isfield(options,'output')
    outp = options.output;
end

switch outp
    case 1
        ds = score;
    case 2
        ds = [ds score];
    case 3
        ds = [ds (Bft_sigmoid(score,options.param)-0.5)*2];
end
