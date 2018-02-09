function [pseudoperm] = LmoPermuTiebreak(c)
% minimize c'*x over permutahedron with tiebreaking constraint

n = length(c);
[cp,perm] = sort(c);
[~,sig] = sort(perm);
z = (n:-1:1)';
pseudoperm=z(sig);
if sig(1) > sig(n) % i.e. z(sig(1)) < z(sig(n))
    return
else
    idxs = [(1:sig(1)-1), (sig(1)+1:sig(n)-1), (sig(n)+1:n)];
    cdrop = cp(idxs);
    d = cdrop - 1/2*(c(sig(1))+c(sig(n)));
    k=find(d<0,1,'last');
    if isempty(k)
        k=1;
    else
        k = k+1;
    end
    noks = [1:k-1,k+2:n];
    znok = z(noks);
    [~,subsig] = sort(perm(idxs));
    pseudoperm = [z(k+1);znok(subsig);z(k)];

end

end
