function [Sproj] = proj2RmatWithLinProg(A)

% remove diagonal
% A = A - diag(diag(A));

n = size(A, 1);
ntri = n*(n+1)/2;
[lti, ltj, ~] = find(tril(ones(n),0));
ltidxs = sub2ind([n, n], lti, ltj);
square2ltmat = sparse((1:ntri), ltidxs, ones(1, ntri), ntri, n^2);

kCons = 1;
iis = []; jjs=[]; vvs=[];

for idiag=0:n-1
    mydiag = (idiag+1) + (n+1)*(0:n-(idiag+1));
    [mydiag,~,~] = find(square2ltmat(:,mydiag));
    mydiag = mydiag';
    n_add = length(mydiag);
    % diagonal entries larger than additional variable underneath
    jjs = [jjs, mydiag, (ntri+idiag+1)*ones(1, n_add)];
    vvs = [vvs, -ones(1, n_add), ones(1, n_add)];
    iis = [iis, kCons:kCons+n_add-1, kCons:kCons+n_add-1];
    kCons = iis(end)+1;
    % diagonal entries smaller than additional variable above
    if idiag>0
        jjs = [jjs, mydiag, (ntri+idiag)*ones(1, n_add)];
        vvs = [vvs, ones(1, n_add), -ones(1, n_add)];
        iis = [iis, kCons:kCons+n_add-1, kCons:kCons+n_add-1];
        kCons = iis(end)+1;
    end
end

nCons = max(iis);
lineqmat = sparse(iis, jjs, vvs, nCons, ntri + n);

slt = square2ltmat*A(:);

% Solve the projection of s with a convex optimization solver
% fprintf('begin call to cvx');
% tic;
% biglineqmatabs1 = [speye(ntri), sparse(ntri,n), -speye(ntri), sparse(ntri,n)];
% biglineqmatabs2 = [-speye(ntri), sparse(ntri,n), -speye(ntri), sparse(ntri,n)];
% % biglineqmatpos = [-speye(ntri+n), sparse(ntri+n,ntri+n)];
% biglineqmatpos = [-speye(ntri+n), -speye(ntri+n)];
% biglineqmatcons = [lineqmat, sparse(nCons, ntri+n)];
% biglineqmat = [biglineqmatabs1; biglineqmatabs2; biglineqmatpos; biglineqmatcons];
% bigf = sparse([zeros(ntri+n,1), ones(ntri+n,1)]);
% bigb = sparse([slt; -slt; zeros(nCons + ntri + n,1)]);
% 
% xsol = linprog(bigf,biglineqmat,bigb);


% Is this faster ?
biglineqmatabs1 = [-speye(ntri), sparse(ntri,n), speye(ntri), sparse(ntri,n)];
biglineqmatabs2 = [-speye(ntri), sparse(ntri,n), -speye(ntri), sparse(ntri,n)];
% biglineqmatpos = [-speye(ntri+n), sparse(ntri+n,ntri+n)];
% biglineqmatpos = [-speye(ntri+n), -speye(ntri+n)];
biglineqmatcons = [sparse(nCons, ntri+n), lineqmat];
biglineqmat = [biglineqmatabs1; biglineqmatabs2; biglineqmatcons];
bigf = sparse([ones(ntri+n,1), zeros(ntri+n,1)]);
bigb = sparse([slt; -slt; zeros(nCons,1)]);
options = optimoptions('linprog','Algorithm','dual-simplex');%,'OptimalityTolerance',1e-7);
xsol = linprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri+n),1),[], options);
% xsol = linprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri+n),1),[]);

sprojvec = square2ltmat'*xsol(ntri+n+1:2*ntri+n);


% sprojvec = square2ltmat'*xsol(1:ntri);
Sproj = reshape(sprojvec, n, n);
Sproj = Sproj + tril(Sproj,-1)';

