function [mys] = projconssparse(b,a,sparam)

p = length(b);
if nargin <= 2 || sparam > p
    sparam = p;
end
[bsrt,myperm] = sort(b, 'descend');
bfirst = bsrt(1:sparam);
myyfirst = projconssorted(bfirst,a);
mys = zeros(1,p);
mys(myperm(1:sparam)) = myyfirst;

end


function [myy] = projconssorted(bb,aa)
% provide a vector b sorted by decreasing value
pp = length(bb);
bsum = sum(bb);

if aa > bsum
    myy = bb + 1/pp * (aa - bsum);
% else
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