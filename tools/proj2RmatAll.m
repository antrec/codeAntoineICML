function [Sproj] = proj2RmatAll(A)

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
elnetpar = 0;
cvx_begin quiet
    variable x(ntri + n);
    minimize elnetpar*norm(x(1:ntri) - slt, 2) + (1-elnetpar)*norm(x(1:ntri) - slt, 1)
    subject to
        lineqmat*x <= 0;
        x >= 0;
cvx_end
% toc;
% Get back to matrix form
sprojvec = square2ltmat'*x(1:ntri);
Sproj = reshape(sprojvec, n, n);
Sproj = Sproj + tril(Sproj,-1)';
% Sproj = Sproj - diag(diag(Sproj));
