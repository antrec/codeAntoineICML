%Calculates the positional distance matrix for the ppc-based fusion:
%K(i,j) is sum of weighted distances between the positions of p(i) and p(j),
%for all permutations p passed as columns of P, with weights w

%JY Goulermas, 2013

function [K] = pdm(P, w)
  [n,b] = size(P);
  P     = dma.invperm(P); %need positions on column vectors
  K     = zeros(n);
  for k = 1 : b
    p = P(:,k);
    A = bsxfun(@minus, p, p'); %o=ones(n,1); A= p*o'-o*p'
    A = A .^ 2;  %as in paper
    %A = abs(A); %an alternative
    K = K + w(k) * A;
  end
end






































