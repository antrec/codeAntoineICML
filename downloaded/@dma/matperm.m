%Converts a permutation from a vector/matrix to a matrix/vector representation.
%No input validation is performed!

%JY Goulermas, 2013

function [out] = matperm(in)

  [n,b] = size(in);

  if isvector(in) %permutation vector given
    out = sparse(1:max(n,b), in, 1);
  else %permutation matrix given
    [out,~] = find(in');
  end

end

