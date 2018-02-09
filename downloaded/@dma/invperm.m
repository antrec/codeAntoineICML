%Generates the inverse permutation vector ip of a given permutation vector p,
%i.e. ip(p)= [1:n]
%If p is an n-length permutation vector, it can be a row or column vector.
%If p is an nxb matrix of b permutation vectors, the operation is applied column-wise.

%JY Goulermas

function [ip] = invperm(p)

  ip = NaN(size(p));

  if isvector(p)
    ip(p) = 1 : length(p);
  else
    for k = 1 : size(p,2)
      ip( p(:,k), k ) = 1 : size(p,1);
    end
  end

end