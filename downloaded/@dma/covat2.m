%From paper: Havens et al., A new formulation of the coVAT algorithm, IJOIS, 2012
%
%X      : an nxm matrix with the two-mode data.
%perm   : the permutation vectors

%JY Goulermas, 2013

function [perm] = covat2(X)

  Dr       = squareform( pdist(X,  'euclidean') );
  Dc       = squareform( pdist(X', 'euclidean') );
  perm.row = dma.vat(Dr);
  perm.col = dma.vat(Dc);
 
end