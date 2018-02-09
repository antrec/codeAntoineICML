%Converts a ranking p (defined as an inverse permutation vector; e.g., Critchlow 84),
%which assigns object with p(i) rank i,
%to a binary Relation matrix R, where R(i,j) = 1 iff p ranks object i before object j.

%A Kostopoulos & JY Goulermas, 2013

function [R] = relation(p)

  n = max(size(p));
  R = zeros(n);
  %assert( isvector(p) && isequal( unique(p(:)), (1:n)') )
  
  for i = 1 : n
    idx          = p(i+1:end);
    R(p(i), idx) = ones(1,n-i);
  end
  
end

%--------------------------------------------------------------------------------------------------

%older version not used
function [R] = relation_naive(p)

  n = max(size(p));
  p = dma.invperm(p);
  R = zeros(n);
  for i=1:n
    for j=1:n
      if p(i) < p(j)
        R(i,j) = 1;
      end
    end
  end

end

%--------------------------------------------------------------------------------------------------