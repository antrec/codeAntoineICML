function [Sproj] = proj2Rmat2(A, opts)

n = size(A,1);
opts_def.doWeakR = false;
opts_def.kMax = n;
opts_def.ubs = inf*ones(1,n);

if nargin == 1
    opts = opts_def;
else
    opts = build_opts(opts_def,opts);
end
kMax    = opts.kMax;
doWeakR = opts.doWeakR;
ubs     = opts.ubs;




n = size(A, 1);
ntri = n*(n+1)/2;
[lti, ltj, ~] = find(tril(ones(n),0));
ltidxs = sub2ind([n, n], lti, ltj);
square2ltmat = sparse((1:ntri), ltidxs, ones(1, ntri), ntri, n^2);

kCons = 1;
iis = []; jjs=[]; vvs=[];

for idiag=0:n-1
    myDiagIdx = (idiag+1) + (n+1)*(0:n-(idiag+1));
    [myDiagIdx,~,~] = find(square2ltmat(:,myDiagIdx));
    myDiagIdx = myDiagIdx';
    n_add = length(myDiagIdx);

    % diagonal entries larger than additional variable underneath
    jjs = [jjs, myDiagIdx, (ntri+idiag+1)*ones(1, n_add)];
    vvs = [vvs, -ones(1, n_add), ones(1, n_add)];
    iis = [iis, kCons:kCons+n_add-1, kCons:kCons+n_add-1];
    kCons = iis(end)+1;

    % diagonal entries smaller than additional variable above
    if idiag>0
        jjs = [jjs, myDiagIdx, (ntri+idiag)*ones(1, n_add)];
        vvs = [vvs, ones(1, n_add), -ones(1, n_add)];
        iis = [iis, kCons:kCons+n_add-1, kCons:kCons+n_add-1];
        kCons = iis(end)+1;
    end
end                                                                                                                                                                                                                                                                                                       

nCons = max(iis);
lineqmat = sparse(iis, jjs, vvs, nCons, ntri + n);

slt = square2ltmat*A(:);



% Is this faster ?
biglineqmatabs1 = [-speye(ntri), speye(ntri), sparse(ntri,n)];
biglineqmatabs2 = [-speye(ntri), -speye(ntri), sparse(ntri,n)];
% biglineqmatpos = [-speye(ntri+n), sparse(ntri+n,ntri+n)];
% biglineqmatpos = [-speye(ntri+n), -speye(ntri+n)];
biglineqmatcons = [sparse(nCons, ntri), lineqmat];
biglineqmat = [biglineqmatabs1; biglineqmatabs2; biglineqmatcons];
bigf = sparse([ones(ntri,1); zeros(ntri+n,1)]);
bigb = sparse([slt; -slt; zeros(nCons,1)]);
options = optimoptions('linprog','Algorithm','dual-simplex');%,'OptimalityTolerance',1e-7);

%     ub = [max(slt(:))*ones(ntri+n,1); ubval*ones(ntri+n,1)];
ub = [max(ubs)*ones(1,2*ntri), ubs];
% xsol = linprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri),1),ub, options);
xsol = callLinprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri),1),ub, options);
%     xsol = linprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri+n),1),[], options);
% xsol = linprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri+n),1),[]);

sprojvec = square2ltmat'*xsol(ntri+1:2*ntri);


% sprojvec = square2ltmat'*xsol(1:ntri);
Sproj = reshape(sprojvec, n, n);
Sproj = Sproj + tril(Sproj,-1)';


    
end
