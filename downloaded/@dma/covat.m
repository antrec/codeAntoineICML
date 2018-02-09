%From paper: Bezdek et al., Visual assessment of cluster tendency for rectangular dissimilarity matrices, TFS, 2007
%
%X      : an nxm matrix with the two-mode data.
%perm   : the permutation vectors
%details: a structure that includes the (ordered) Dr, Dc & DrUc matrices and original permutation vector from VAT

%JY Goulermas, 2010

function [perm, details] = covat(X)

  [m,n] = size(X);

  %initialise
  Dr    = squareform( pdist(X,  'euclidean') );
  Dc    = squareform( pdist(X', 'euclidean') );
  d_avg = mean(X(:)); 
  sc_r  = d_avg / ( sum(Dr(:)) / (m^2-m) ); %scale Dr & Dc to adjust contrast
  sc_c  = d_avg / ( sum(Dc(:)) / (n^2-n) );
  Dr    = Dr * sc_r;
  Dc    = Dc * sc_c;
  Drc   = [Dr, X; X', Dc];  %Drc=[inf(size(Dr)), X; X', inf(size(Dc))]; %for testing (just model the bipartite)
  P     = dma.vat(Drc);

  perm.row    = P(P<=m);
  perm.col    = P(P> m) - m;
  details.Drc = Drc(P,P);
  details.Dr  = Dr (perm.row,perm.row); %approximations of VAT(Dr) & VAT(Dc)
  details.Dc  = Dc (perm.col,perm.col);
  
end