function [Sproj] = proj2RmatSparse(A, k, square2ltmat, lineqmat)
% Projection on R matrices

n = size(A, 1);
if (nargin == 1 || nargin == 3)
    k=n;
end
if nargin <= 2
    [square2ltmat, lineqmat] = buildRConsSparse(n, k);
end

nel = (k+1)*(2*n-k)/2;
A = tril(A,k) - tril(A,-k-1);
slt = square2ltmat*A(:);

% Solve the projection of s with a convex optimization solver
tic;
fprintf('begin call to cvx');
cvx_begin quiet
    variable x(nel + k+1);
    
    minimize norm(x(1:nel) - slt, 1)
    subject to
        lineqmat*x <= 0;
        x >= 0;
cvx_end
toc;

% Get back to matrix form
sprojvec = square2ltmat'*x(1:nel);
Sproj = reshape(sprojvec, n, n);
Sproj = Sproj + tril(Sproj,-1)';
% Sproj = Sproj - diag(diag(Sproj));
