function [Sproj] = proj2RmatWithLinProg(A,doWeakR, ubval)

% remove diagonal
% A = A - diag(diag(A));
if nargin > 1 && doWeakR
        
        % generate constraints for R matrices
        n = size(A,1);
        
        
        
        
        
        
        
nCons = n*(n+1);
ntri = n*(n+1)/2;
iis = [];
jjs = [];
vvs = [];

[idiag,jdiag,~]=find(eye(n));
diagIdx = sub2ind([n,n],idiag,jdiag);
idxConsMax = 0;

for iCol=1:n-1
    % columnwise R-constraints ( A_{i+1,j} <= A{i,j} in lower triangle)
    iIdxUp = diagIdx(iCol):diagIdx(iCol)+(n-iCol-1);
    iIdxDown = diagIdx(iCol)+1:diagIdx(iCol)+(n-iCol);
    idxCons = idxConsMax + (1:n-iCol);
    jjs = [jjs, idxCons, idxCons];
    iis = [iis, iIdxUp, iIdxDown];
    vvs = [vvs, -ones(1,n-iCol), ones(1,n-iCol)];
    idxConsMax = idxConsMax + n-iCol;
    
    %linewise R-constraints ( A_{i,j} <= A_{i,j+1} in lower triangle)
    iIdxLeft = diagIdx(iCol) + (1:(n-iCol));
    iIdxRight = diagIdx(iCol+1) + (0:(n-iCol-1));
    idxCons = idxConsMax + (1:n-iCol);
    jjs = [jjs, idxCons, idxCons];
    iis = [iis, iIdxLeft, iIdxRight];
    vvs = [vvs, ones(1,n-iCol), -ones(1,n-iCol)];
    idxConsMax = idxConsMax + n-iCol;

end

[lti, ltj, ~] = find(tril(ones(n),0));
ltidxs = sub2ind([n, n], lti, ltj);
square2ltmat = sparse((1:ntri), ltidxs, ones(1, ntri), ntri, n^2);
        
lineqmat = sparse(jjs,iis,vvs,nCons,n^2);
lineqmat = lineqmat*square2ltmat';
        
%         
%         
%         
%         
%         
%         nCons = (n-1)*(n-2);
%         ntri = n*(n-1)/2;
%         tmpvec = [1,floor(2:0.5:n-1)];
%         iis = [];
%         jjs = [];
%         vvs = [];
%         itmpval = 0;
%         jtmpval = 0;
% 
%         csvec = cumsum([0,2:n]);
%         jproj=[];
%         jprec=0;
% 
%         for icol=1:n-2
% 
%             % constraints of type R_{i,j} <= R_{i,j+1}
%             icoords = itmpval + tmpvec(1:2*(n-1-icol));
%             iis = [iis, icoords];
%             itmpval = itmpval + (n-icol);
%             jcoords = jtmpval + floor(1:0.5:(n-icol-0.5));
%             jjs = [jjs, jcoords];
%             jtmpval = jtmpval + (n-1-icol);
%             altsign = 2*(mod((1:2*(n-1-icol)),2)-1/2);
%             vvs = [vvs, altsign];
% 
%             % also constraints of type R_{i+1,j} <= R_{i,j}
%             thaticoords = (icol+1:n-1) + (icol-1)*n - csvec(icol);
%             jcoords = jtmpval + (1:n-icol-1);
%             vvs = [vvs, -ones(1, n-icol-1)];
%             jjs = [jjs, jcoords];
%             iis = [iis, thaticoords];
% 
%             thaticoords = (icol+1:n-1) + (icol)*n - csvec(icol+1);
%             vvs = [vvs, ones(1, n-icol-1)];
%             jjs = [jjs, jcoords];
%             iis = [iis, thaticoords];
% 
%             jtmpval = jtmpval + (n-1-icol);
% 
%             % also construct projection from vectorized matrix n^2 to only the
%             % entries in the lower triangle (n*(n-1)/2))
%             jproj = [jproj, jprec+icol+(1:n-icol)];
%             jprec = jproj(end);
% 
%         end
% 
%         % finish constructing the projection from vectorized matrix n^2 to only the
%         % entries in the lower triangle (n*(n-1)/2))
%         icol=n-1;
%         jproj = [jproj, jprec+icol+(1:n-icol)];
%         vproj = ones(1, length(jproj));
%         iproj = (1:length(jproj));
%         square2ltmat = sparse(iproj, jproj, vproj, ntri, n^2);
%         
%         % Constraint matrix. Note that I used iis for columns and jjs for rows by 
%         % mistake...
%         lineqmat = sparse(jjs, iis, -vvs, nCons, ntri);
        
        
        % Vectorize the input S and keep only entries in the lower triangle
        slt = square2ltmat*A(:);

        
        biglineqmatabs1 = [-speye(ntri), speye(ntri)];
        biglineqmatabs2 = [-speye(ntri), -speye(ntri)];
        biglineqmatcons = [sparse(nCons, ntri), lineqmat];
        
%         % Test heuristic : add upper bound on diagonal
%         biglineq
        
        biglineqmat = [biglineqmatabs1; biglineqmatabs2; biglineqmatcons];
        bigf = sparse([ones(ntri,1), zeros(ntri,1)]);
        bigb = sparse([slt; -slt; zeros(nCons,1)]);
        
        % Call cvx opt. solver
        options = optimoptions('linprog','Algorithm','dual-simplex');%,'OptimalityTolerance',1e-7);
        % Add upper bound ?
        ub = [slt; ubval*ones(ntri,1)];
%         ub = ubval*ones(2*ntri,1);
%         xsol = linprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri),1),ub, options);
        xsol = callLinprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri),1),ub, options);
%         xsol = linprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri),1),[], options);

        % Get back to matrix form
        sprojvec = square2ltmat'*xsol(ntri+1:end);
        Sproj = reshape(sprojvec, n, n);
        Sproj = Sproj + Sproj' - diag(diag(Sproj));
        
        return;
            
    
else

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
    
    ub = [max(slt(:))*ones(ntri+n,1); ubval*ones(ntri+n,1)];
%     xsol = linprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri),1),ub, options);
    xsol = callLinprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri),1),ub, options);

%     xsol = linprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri+n),1),[], options);
    % xsol = linprog(bigf,biglineqmat,bigb, [], [], sparse(2*(ntri+n),1),[]);

    sprojvec = square2ltmat'*xsol(ntri+n+1:2*ntri+n);


    % sprojvec = square2ltmat'*xsol(1:ntri);
    Sproj = reshape(sprojvec, n, n);
    Sproj = Sproj + tril(Sproj,-1)';

    
    
end
