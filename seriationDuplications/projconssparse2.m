function [mys] = projconssparse2(b,a, opts)


p = length(b);
opts_def = [];
opts_def.sparam = p;
opts_def.ub = inf;
if nargin <= 2
    opts = opts_def;
else
    opts = build_opts(opts_def,opts);
end
sparam = opts.sparam;
sparam = min(p,sparam);
ub = opts.ub;


[bsrt,myperm] = sort(b, 'descend');
bfirst = bsrt(1:sparam);
% if sparam <= 1
%     fprintf('%d',sparam);
% end
myyfirst = projconssorted(bfirst,a, ub);
mys = zeros(1,p);
mys(myperm(1:sparam)) = myyfirst;

end


function [myy] = projconssorted(bb,aa, ub)
% provide a vector b sorted by decreasing value
pp = length(bb);
bsum = sum(bb);

if aa > bsum
    myy = bb + 1/pp * (aa - bsum);
    if myy(1) > ub
        firstIn = find(myy<ub,1,'first');
        if isempty(firstIn)
            myy = ub*ones(size(myy));
        else
            myy(1:firstIn-1) = ub;
            myy(firstIn:end) = projconssorted(bb(firstIn:end),aa-(firstIn-1)*ub,ub);
        end
    end
else
%     [sbvals,thisp] = sort(b, 'descend');
    bma = bb - (1./(1:pp)') .* (cumsum(bb) - aa);
    fneg = find(bma<0,1,'first');
    if isempty(fneg)
        fneg = length(bma);
    else
        fneg = fneg - 1;
    end
    myvals = zeros(1,pp);
    myvals((1:fneg)) = bb(1:fneg) - 1/fneg * (sum(bb(1:fneg)) - aa);
    myy = myvals;

end

end