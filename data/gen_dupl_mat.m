function [A, Z, c] = gen_dupl_mat(S, sz_ratio, dupl_prop)
% generate matrix with duplications from matrix S. If S is pre-R, this is
% like JPV problem.
if nargin <= 2
    dupl_prop = 0.05;
end
N = length(S);
% dupl_max = 3;
n = floor(N*sz_ratio);
n_dupl = floor(dupl_prop*n);
dupl_idx = randperm(n, n_dupl);
Z = sparse(1:n, 1:n, ones(1,n), N, n);
rep_idx = dupl_idx(randi(n_dupl, 1, N-n));
Z = Z + sparse(n+1:N, rep_idx, ones(1, N-n), N, n);

rpN = randperm(N)';
rpDupl = [(1:n), randi(n,1,N-n)];
rpDupl = rpDupl(randperm(N)');
% rpDupl = randi(n,1,N);
Z = sparse(rpDupl, rpN, ones(1,N), n, N);
Z = Z';

    



A = Z'*S*Z;
c = Z'*ones(N, 1);

% 
% if nargin == 1
%     c = randi(dupl_max-1, n, 1) +1;
%     r = rand(1, n, 1);
%     c(r>dupl_prop) = 1;
% end
