function [square2ltmat, lineqmat] = buildRConsSparse(n, k)
% build inequality constraint matrix for projection onto R matrices
if nargin <= 1
    k=n;
end

nel = (k+1)*(2*n-k)/2;
T = tril(ones(n), -1) - tril(ones(n), k);
[lti, ltj, ~] = find(T');
ltidxs = sub2ind([n, n], lti, ltj);
square2ltmat = sparse((1:nel), ltidxs, ones(1, nel), nel, n^2);

kCons = 1;
iis = []; jjs=[]; vvs=[];

for idiag=0:k-1
    mydiag = (idiag+1) + (n+1)*(0:n-(idiag+1));
    [mydiag,~,~] = find(square2ltmat(:,mydiag));
    mydiag = mydiag';
    n_add = length(mydiag);
    % diagonal entries larger than additional variable underneath
    jjs = [jjs, mydiag, (nel+idiag+1)*ones(1, n_add)];
    vvs = [vvs, -ones(1, n_add), ones(1, n_add)];
    iis = [iis, kCons:kCons+n_add-1, kCons:kCons+n_add-1];
    kCons = iis(end)+1;
    % diagonal entries smaller than additional variable above
    if idiag>0
        jjs = [jjs, mydiag, (nel+idiag)*ones(1, n_add)];
        vvs = [vvs, ones(1, n_add), -ones(1, n_add)];
        iis = [iis, kCons:kCons+n_add-1, kCons:kCons+n_add-1];
        kCons = iis(end)+1;
    end
end

nCons = max(iis);
lineqmat = sparse(iis, jjs, vvs, nCons, nel + k+1);