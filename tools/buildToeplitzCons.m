function [lineqmat, eqmat, square2ltmat] = buildToeplitzCons(n)


% % constraints matrix
nCons = n*(n+1)/2;
ntri = n*(n+1)/2;
% 
% lineqmat = sparse(jcons, icons, vcons, nCons, n*(n-1)/2);
% % lineqmat = sparse(jcons, icons, vcons, nCons, n^2);
% % lineqmat = lineqmat*square2ltmat';

% equality constraints
kCons = 1;
iis = []; jjs=[]; vvs=[];
nextdiagidx = 1 + (n+1)*(0:n-1);
ltidxs = zeros(1, ntri);
ltidxs(1:n) = nextdiagidx;
endidx = n;
for idiag=1:n-1
    mydiagidx = nextdiagidx;
    nextdiagidx = (idiag+1) + (n+1)*(0:n-(idiag+1));
    jjs = [jjs, mydiagidx(1:end-1), mydiagidx(2:end)];
    n_add = length(mydiagidx)-1;
    vvs = [vvs, -ones(1, n_add), ones(1, n_add)];
    iis = [iis, kCons:kCons+n_add-1, kCons:kCons+n_add-1];
    kCons = iis(end)+1;
    ltidxs(endidx+1:endidx+length(nextdiagidx)) = nextdiagidx;
    endidx = endidx+length(nextdiagidx);
end
Cbig = sparse(iis, jjs, vvs, nCons, n^2);
square2ltmat = sparse((1:ntri), ltidxs, ones(1, ntri), ntri, n^2);
eqmat = Cbig*square2ltmat';

% inequality constraints
firstcol = (1:n);
jjs = [firstcol(1:end-1), firstcol(2:end)];
vvs = [-ones(1, n-1), ones(1, n-1)];
iis = [(1:n-1), (1:n-1)];
Cbig = sparse(iis, jjs, vvs, n-1, n^2);
lineqmat = Cbig*square2ltmat';
