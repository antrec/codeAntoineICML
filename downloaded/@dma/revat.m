%From paper: Huband et al., reVAT, revised visual assessment of cluster tendency, NAFIPS, 2004
%
%D      : a nxn dissimilarity matrix
%m      : the sample size for bigvat
%r      : the fraction of the row average used for thresholding (default=0.5)
%perm   : a permutation vector
%pivot  : the pattern indices chosen to form neighbourhoods
%profile: used for the scaled (highest=1) bar plots 

%JY Goulermas, 2009

function [perm, details] = revat(D, m, r)

  if ~exist('r', 'var')
    r = 0.5;
  end

  n    = size(D,1);
  perm = [];
  out  = true(1,n) ;
  iter = 0;
  
  while sum(out)
    iter                = iter + 1;
    i                   = find(out, 1, 'first');
    %idx = find(out); i = floor( random('unif', 1, numel(idx), 1) ); i = idx(i); %random
    pivot(iter)         = i;
    row                 = D(i,:);
    delta               = mean(row) * r;
    tall_set{iter}      = find( row < delta & out ); %thresholding
    perm                = [perm; tall_set{iter}'];
    out(tall_set{iter}) = 0;
  end
  
  details.pivot   = pivot;
  details.profile = D(pivot, perm);
  details.profile = 1 - details.profile ./ max(D(:)); 
  perm            = BigvatIndices(m);
    
  
  %----------------------------------------------
  %From paper: Huband et al., bigVAT: visual assessment of cluster tendency for large datasets, Pattern Recogition, 2005
  function [p] = BigvatIndices(new_size)
   
    if new_size <= 1
      new_size = new_size * n;
    end

    if n <= new_size
      p = perm;
    else
      p = [];
      for i = 1 : numel(tall_set)
        portion = floor( numel(tall_set{i}) * new_size / n );
        shuffled = shuffle( tall_set{i} );
        selected = shuffled(1:portion);
        p        = [p; selected'];
      end
    end

  end
  %----------------------------------------------
  
end