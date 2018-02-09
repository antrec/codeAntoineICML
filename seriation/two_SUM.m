function [f, df, ddf] = two_SUM(x, A)

n = length(A);
% c = 1./2*(n+1)*ones(n,1); % center of the convex hull of permutations.
doTrans=false;
if size(x) ~= [n,1]
    if size(x) == [1,n]
        x = x';
        doTrans = true;
    else
        fprintf('x must be of length n !');
    end
end
        

% Compute Laplacian
if issparse(A)
    L = spdiags(sum(A,2), 0, n, n) - A;
    L= L+(1e-10)*speye(n);
else
    L = diag(sum(A,2)) - A;
    L = L+(1e-10)*eye(n);
end


f = x'*L*x;
df = 2*L*x;
if doTrans
    df = df';
end
ddf = L;

if issparse(A)
    f = full(f);
    df = full(df);
end

end