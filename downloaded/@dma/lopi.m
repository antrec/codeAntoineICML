%Locally Optimal Pairwise Interchange (LOPI), from paper: Brusco & Stahl 2000
%It minimises single- or multi-commodity QAPs expressed as: Sum_i trace( P* D_i *P' * F_i ).
%perm   : the optimising permutation
%details: history information

%JY Goulermas, 2014

function [perm, details] = lopi(D, F, p)
  if ~iscell(D)
    D = {D};
    F = {F};
  end
  time = cputime;
  t    = 0;
  stop = false;
  while ~stop
    t = t + 1;
    [best(t), i, j] = LocalOptimisation(D, F, p);
    if best(t) < 0 %further improvement possible
      temp = p(i); p(i) = p(j); p(j) = temp;
      fprintf('\n LOPI iteration %d, with best delta %6.4f, swapping indices (%d,%d)', t, best(t), i, j);
    else
      stop = true;
      cost  = 0;
      for k = 1 : numel(D)
        cost = cost + sum(sum( D{k}(p,p) .* F{k} ));
      end
      fprintf('\n LOPI iteration %d, with best delta %6.4f, and final cost %6.4f', t, best(t), cost);
    end
  end
  perm = p;
  details.time  = cputime - time;
  details.delta = best;
  details.cost  = cost;
  details.evals = t * nchoosek(size(D{1},1),2);
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Heuristic search for finding the best (i,j) swap. If best_delta<0, then an
%improving swap has been found, otherwise one does not exist
function [best_delta, best_i, best_j] = LocalOptimisation(D, F, p)
  n          = size(D{1},1);
  L          = numel(D);
  best_delta = inf;
  for i = 2 : n
    for j = 1 : i-1
      idx        = [1:i-1, i+1:j-1, j+1:n]; %indices affected by the swap
      delta_swap = @(D,F) 2*( -D(p(i),p(idx)) +D(p(j),p(idx)) ) * ( F(i,idx) - F(j,idx) )'; %=new-old
      delta = 0;
      for k = 1 : L
        delta = delta + delta_swap(D{k},F{k});
      end
      if delta < best_delta %better than best so far
        best_delta = delta;
        best_i     = i;
        best_j     = j;
      end
    end
  end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~














