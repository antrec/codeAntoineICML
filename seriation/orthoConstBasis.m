function [B] = orthoConstBasis(n)
% basis of the hyperplane orthogonal to constant vector in R^n

% construct orthogonal basis
I = eye(n);
B = cumsum(I, 2);
D = diag(0:n-1);
B(:,1:n-1) = B(:,1:n-1) - D(:, 2:n);

% normalize columns
vec = (1:n);
vec(1:n-1) = vec(1:n-1).*vec(2:n);
N = diag(1./sqrt(vec));
B = B*N;