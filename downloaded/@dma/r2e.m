%From paper: Chen, Generalized Association Plots, Stat. Sinica, 2002.
%
%Applies the rank-two ellipse (R2E) seriation method which is based on the iterative application
%of the correlation function onto the samples stored as rows in D.
%X      : an nxm matrix
%mode   : 1 or 2 to specify the operating mode
%perm   : the permutation vector(s)
%details: a structure to return the number of iterations & projections

%JY Goulermas, 2009

function [perm, details] = r2e(X, mode)

  %Repeat the correlation procedure
  D     = X;
  iters = 0;
  while rank(D) > 2 && iters < 2000
    iters = iters + 1;
    D     = corr(D', 'type', 'pearson'); %faster than: D = 1 - squareform( pdist(D, 'correlation') );
    %note(2017): corr(of a symmetric) is not always exactly symmetric, and this causes problems with CONCOR's fixed points or rank>1
  end%while                                         

  %Uses the projections Y of all points D onto the first two eigenvectors as D*[m1,m2] = [m1,m2]*diag([l1,l2])
  %which lie on a 2D ellipse with kernel inv(L)
  [M,L]    = eig(D);
  [~, idx] = sort(diag(L), 'descend');
  idx      = idx([1,2]);
  Y        = M(:,idx) * L(idx,idx); %scaling with L is not really needed as the order of angles cannot change
  
  %Order points by tracking the angles of the points
  angles = atan2( Y(:,2), Y(:,1) );
  perm   = CircularSort(angles);

  details.iters = iters;
  details.Y     = Y;
  
  if mode == 2
    [perm2,details2] = dma.r2e(X',1);
    perm             = struct('row', perm, ...
                              'col', perm2 ...
                             );
    details          = struct('iters_row', details.iters, ...
                              'iters_col', details2.iters, ...
                              'Y_row',     details.Y, ...
                              'Y_col',     details2.Y ...
                             );
  end%mode

end
