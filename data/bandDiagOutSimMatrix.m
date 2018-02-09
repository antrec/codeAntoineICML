function [simMatrix] = bandDiagOutSimMatrix(n, kBD, nOut)
% generate similarity matrix S = D + N where D is band diagonal of size 
% kBD and N is nOut sparse out-of-band-diagonal terms

simMatrix = tril(ones(n) - tril(ones(n),-kBD));% - triu(ones(n),+kBD);
simMatrix = sparse(simMatrix);
% potential indices for putting nOut
[iin,jin,~] = find(tril(~simMatrix));
idxs = sub2ind([n n], iin, jin);
nOutHalf = floor(nOut); % divide nOut by 2 to fill lower triangle
nidxs = length(idxs);
% Chose nOutHalf indexes among the nidxs available
if nOutHalf > nidxs
    fprintf('nOut too large : matrix full of ones !');
    simMatrix = ones(n);
    return;
end
thisIdxs = randperm(nidxs, nOutHalf);
thisIdxs = idxs(thisIdxs);
[iout, jout] = ind2sub([n n], thisIdxs);
simMatrix = simMatrix + sparse(iout, jout, ones(size(iout)), n, n);

% Symmetrize
simMatrix = simMatrix + tril(simMatrix,-1)';

end