function [f, df] = transfo2permhyperplane(U, y, f_handle)

[n, m] = size(U);
if m ~= n-1
    fprintf('U should be of size n x n-1 !');
end
if length(y) ~= n-1
    fprintf('y should be of length n-1 !');
end

% to make things easier let's add the constant to be in same hyperplane as
% the permutahedron
c = 1./2*(n+1)*ones(n,1); % center of the convex hull of permutation vectors.
x = U*y + c;
[f, df] = f_handle(x);
df = U'*df;