function [Xtest] = proj2DuplConsSparseClosedForm(S, Z, A, opts)
% projects matrix S onto set of constraints of J.P. V. Seriation with Duplications
% problem. Z : nxN assignment matrix, S : NxN to be R mat, A : nxn dupl.
% mat.

n = size(S,1);
smalln = size(A,1);

% keep only lower triangular part for symmetric matrix
if issymmetric(A)
    A = tril(A,0);
elseif istriu(A)
    A = A';
elseif ~istril(A)
    fprintf('A is not symmetric nor triangular !');
end

if issymmetric(S)
    S = tril(S,0);
elseif istriu(A)
    S = S';
elseif ~istril(S)
    fprintf('S is not symmetric nor triangular !');
%     imagesc(S); colorbar;
%     sum(sum(abs(tril(S) - triu(S)')))
%     sum(sum(abs(tril(S))))
%     pause;
%     imagesc(tril(S)+tril(S,-1)' - triu(S) - triu(S,1)'); colorbar;
%     pause;
end

[is,js,vs]=find(A);
nbnz = length(is);
difj = js(2:end)-js(1:end-1);
switchj=find(difj);
switchj = switchj+1;
switchj = [switchj; length(js)+1];

opts_def.batchSize = nbnz;
opts_def.softCons = false;
opts_def.sparsprop = 1;
opts_def.elnet = 0.5;
opts_def.dh = n;
opts_def.ub = inf;
% opts_def.maxval = inf;
if nargin == 3
    opts = opts_def;
else
    opts = build_opts(opts_def, opts);
end
batchSize = opts.batchSize;
softCons = opts.softCons;
elnet = opts.elnet;
sparsprop = opts.sparsprop;
dh = opts.dh;
ubval = opts.ub;

Svec = S(:);

if softCons
    newvals = Svec;
else
    newvals = zeros(1,n^2);
end
% 
% conscolsub = [];
% consrowsub = [];

nextswitchidx = 1;
nextswitch = switchj(nextswitchidx);
jj = js(1);
ls = find(Z(jj, :));
for nzi=1:nbnz
    
    ii = is(nzi);
    ks = find(Z(ii, :));
    if nzi == nextswitch
        nextswitchidx = nextswitchidx + 1;
        nextswitch = switchj(nextswitchidx);
        jj = js(nzi);
        ls = find(Z(jj, :));
    end
    [K,L] = meshgrid(ks, ls);
    kk = K(:);
    ll = L(:);
    % Keep only inds from lower triangle
    kkk = max(kk, ll);
    lll = min(kk, ll);
    myinds = sub2ind([n n], kkk, lll);
%     % OR KEEP ONLY VALUES IN A GIVEN BAND ?!
%     inband = find(kkk-lll <= dh);
%     if isempty(inband)
%         [~,inband]=min(kkk-lll);
%     end
%     myinds = sub2ind([n n], kkk(inband), lll(inband));
    
    % Compute new values projected on constraints set
    ninds = length(myinds);
    bvals = Svec(myinds);
    aij = vs(nzi);

    if ninds == 1
        newvals(myinds) = aij;%1e3;%bvals;
%         continue
    else
        % Test something : have a target band-size for S
        sparam = ninds;
%         sparam = max(1,sum(kkk-lll <= dh));
%         sparam = min(ninds,max(1,floor(ninds*sparsprop)));
%         sparam = min(sparam, nnz(bvals));
%         sparam = max(2,sparam);
%         sparam = max(sparam, max(length(ks),length(ls)));
%         sparam = ninds;
%         newvals(myinds) = projconssparse(bvals, aij, sparam);
        thisopts = [];
        thisopts.sparam = sparam;
        thisopts.ub = ubval;
        newvals(myinds) = projconssparse2(bvals,aij,thisopts);
%         bsum = sum(bvals);
%         if aij >= bsum
%             newvals(myinds) = bvals + 1/ninds * (aij - bsum);
% 
%         else
%             [sbvals,thisp] = sort(bvals, 'descend');
%             bma = sbvals - (1./(1:ninds)') .* (cumsum(sbvals) - aij);
%             fneg = find(bma<0,1,'first');
%             if isempty(fneg)
%                 fneg = length(bma);
%             else
%                 fneg = fneg - 1;
%             end
%             myvals = zeros(1,ninds);
%             myvals(thisp(1:fneg)) = sbvals(1:fneg) - 1/fneg * (sum(sbvals(1:fneg)) - aij);
%     %         myvals = sbvals(1:fneg) - 1/fneg * (sum(sbvals(1:fneg)) - aij);
%     %         [~,invthisp] = sort(thisp);
%             newvals(myinds) = myvals;
%     %         newvals(myinds) = 300;
% 
%         end
    
    end

    
    
%     % add these indices to the list of non-zero indices
%     
%     % add line to constraint matrix
%     conscolsub = [conscolsub, myinds'];
%     consrowsub = [consrowsub, nzi*ones(1, length(myinds))];
    
end

Xtest = reshape(newvals, n, n);
Xtest = Xtest + tril(Xtest,-1)';

% 
% C = sparse(consrowsub, conscolsub, ones(1, length(consrowsub)), nbnz, n^2);
% 
% 
% % keep only non-zero columns (should at least exclude columns corresponding
% % to upper diagonal entries since we restrict A and S to lower triangle 
% if softCons
%     onestri = tril(ones(n),0);
%     onestri = onestri(:);
%     nzc = find(onestri);
% else
%     nzc = find(sum(C,1));
% end
% C = C(:,nzc);
% nnzc = length(nzc);
% 
% svals = S(:);
% 
% xvec = zeros(n^2,1);
% 
% if nargin == 4 && batchSize < nbnz
%     
%     % Make mini batches for CVX
%     
%     myBatch = 1:batchSize:nbnz;
%     if myBatch(end) < nbnz
%         myBatch = [myBatch, nbnz];
%     end
%     nBatch = length(myBatch) -1;
%     thisBatchSize=batchSize;
% 
%     for iBatch=1:nBatch % A.R. : could do in parallel
%         iinBatch = myBatch(iBatch):myBatch(iBatch+1);
%         if iBatch==nBatch
%             thisBatchSize = length(iinBatch);
%         end
% %         Csub = sparse(consrowsub(iinBatch), conscolsub(iinBatch), ones(1, thisBatchSize), thisBatchSize, n^2);        
%         Csub = C(iinBatch,:);
%         nzcSub = find(sum(Csub,1));
%         Csub = Csub(:,nzcSub);
%         nzcSub = nzc(nzcSub);
%         sSub = svals(nzcSub);
%         nnzcSub = length(nzcSub);
%         
%         % Run convex programming software to perform projection
%         cvx_begin
%             variable xSub(nnzcSub);
%             minimize elnet*norm(xSub - sSub, 1) + (1-elnet)*norm(xSub - sSub,2)
%             subject to
%                 Csub*xSub == vs(iinBatch);
%                 xSub >= 0;
%         cvx_end
% 
%         xvec(nzcSub) = xSub;
%         
%     end
% 
% else
%     
%     % Make only one big CVX program with all values where S>0
%     
%     % use only values where S > 0
%     svals = svals(nzc);
%     
%     % Run convex programming software to perform projection
%     cvx_begin
%         variable xvals(nnzc);
%         minimize elnet*norm(xvals - svals, 1) + (1-elnet)*norm(xvals - svals,2)
%         subject to
%             C*xvals == vs;
%             xvals >= 0;
%     cvx_end
% 
%     xvec(nzc) = xvals;
% %     xvec = square2ltmat'*xvals;
% 
%     
% end
% 
% 
% % return to matrix form
% X = reshape(xvec, n, n);
% X = X + tril(X,-1)';
%         

end