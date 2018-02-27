function [d, Rproj] = dist2R(A, normarg)

Rproj = proj2RmatAll(A);
if nargin < 2
    d = norm(Rproj-A, 'fro');
else
    d = norm(Rproj-A, normarg);
end

