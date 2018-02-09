%From paper: M Friendly, Corrgrams, The American Statistician, 2002
%
%D   : an nxn distance matrix. The original algorithm requires a correlation matrix.
%      In such case, D should be its negated version as the two lowest eigenvectors
%      are chosen here (since we process dissimilarity matrices)
%perm: the permutation vector

%JY Goulermas, 2012

function [perm] = corder(D)
  
  [M,L]    = eig(D);
  [~, idx] = sort(diag(L), 'ascend'); %smallest eigenvalues as D is a distance matrix
  x        = M(:, idx(1));
  y        = M(:, idx(2));
  perm     = CircularSort( atan2(y,x) );

  %wrong algo from paper below (quadrants mixed because sign(y) is ignored); atan2 should be used instead as above
  %a       = atan(y./x);
  %a(x<=0) = a(x<=0) + pi;
end




  