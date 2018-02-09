%From paper: Michael J Brusco & Stephanie Stahl, Using Quadratic Assignment Methods to Generate Initial
%            Permutations for Least-Squares Unidimensional Scaling of Symmetric Proximity Matrices
%            Journal of Classification 17:197-223, 2000.
%Code from: Michael J Brusco, Hans-Friedrich Kohn & Stephanie Stahl, Heuristic Implementation of
%           Dynamic Programming for Matrix Permutation Problems in Combinatorial Data Analysis
%           Psychometrika 17(3):503-522, 2008. Code largely based on the original Fortran version.
%Using simulated annealing to optimise equivalent least squares seriation
%in QAP form:
%(QAP w_ij=|i-j| template)
%
%D    : an nxn distance matrix.
%perm : the permutation vector

%A Kostopoulos & JY Goulermas, 2014

function [perm] = arsa(D, varargin)

  if mod(length(varargin),2)
    error('dma:arsa', 'missing parameters')
  end

  n     = size(D,1);
  %default values
  REPS  = 1;     %number of starting random permutations
  TMIN  = 0.1;   %minimum temperature
  COOL  = 0.5;   %cooling rate
  ILOOP = 100*n;

  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'nreps'
        REPS = varargin{i+1};
     case 'tmin'
        TMIN = varargin{i+1};
     case 'cool'
        COOL = varargin{i+1};
     otherwise
        error('dma:arsa', 'unknown parameter: %s', varargin{i})
    end
  end
  
  P     = zeros(REPS,n);
  [c,r] = meshgrid(1:n);
  W     = abs((r-c));     %weight matrix w(i,j) = abs(i-j)
  ZMAX  = sum(sum( W .* D )) / 2;
  perm  = (1:n)';
  
  for i=1:REPS %Calculate REPS number of random permutations
    P(i,:) = randperm(n);
  end
  
  for reps = 1:REPS
    p              = P(reps,:);
    p_best         = p;
    criterion      = sum(sum( W .* D(p,p) )) / 2;
    criterion_best = criterion;
    TMAX           = 1;
    
    %Calculate the maximum temperature for simulation
    for l=1:5000
      [i,j] = random_nodes(1,n,1);
      [delta] = calculateDELTA_swap(W,D,p,i,j);
      if delta < 0
        if abs(delta) > TMAX
          TMAX = abs(delta);
        end
      end
    end
    
    schedule    = floor( (log(TMIN)-log(TMAX))/log(COOL) );
    temperature = TMAX;
    
    for sim=1:schedule
      fprintf('\n#Repetition: %d, Iteration: %d / %d', reps, sim, schedule);
      for k=1:ILOOP
        if rand < 0.5%Randomly swap nodes
          [i,j]   = random_nodes(1,n,1);%pick two random nodes, ordered
          [delta] = calculateDELTA_swap(W,D,p,i,j);
          
          if delta > -eps
            criterion    = criterion + delta;
            [p]          = swap(p,i,j);
            if criterion > criterion_best
              criterion_best = criterion;
              p_best         = p;
            end
          else
            if rand < exp(delta/temperature)
              criterion = criterion + delta;
              [p] = swap(p,i,j);
            end
          end
          
        else%Randomly shift nodes
          [i,j]   = random_nodes(2,n,0);%pick two random nodes, unordered
          [delta] = calculateDELTA_shift(D,p,i,j);
          
          if delta > -eps
            criterion = criterion + delta;
            [p]       = shift(p,i,j);
            if criterion > criterion_best
              criterion_best = criterion;
              p_best         = p;
            end
          else
            if rand < exp(delta/temperature)
              criterion = criterion + delta;
              [p] = shift(p,i,j);
            end
          end
          
        end%if
      end%k
      temperature = temperature * COOL;
    end%sim
    
    if criterion_best > ZMAX
      ZMAX = criterion_best;
      perm = p_best(:);
    end
  end
end

%--------------------------------------------------------------------------------------------------

function [p] = swap(p,i,j)
  %Swap two nodes in a permutation vector
  temp = p(i);
  p(i) = p(j);
  p(j) = temp;
end

function [delta] = calculateDELTA_swap(W,D,p,i,j)
  %Calculate criterion difference when swapping two random nodes
  n     = size(D,1);
  s     = [1:i-1 i+1:j-1 j+1:n];
  delta = ( W(s,i)-W(s,j) )' * ( D(p(s),p(j)) - D(p(s),p(i)) );
end

function [p] = shift(p,i,j)
  %Shift left/right a set of nodes in a permutation vector
  K = p(i);
  if j>i
    p(i:j-1) = p(i+1:j);
  else
    p(i:-1:j+1) = p(i-1:-1:j);
  end
  p(j) = K;
end

function [delta] = calculateDELTA_shift(D,p,i,j)
  %Calculate criterion difference when shifting between two random nodes
  
  if true %if true, uses corrected calculation of delta, if false, uses original fortran calculation of delta
    n    = size(D,1);
    span = abs(j-i);

    if j>i  
      BEG = (1:i-1)';
      MID = (i+1:j)';
      END = (j+1:n)';

      delta = sum(sum( D(p(END),p(MID)) )) - sum(sum( D(p(BEG),p(MID)) ));

      delta = delta + span * sum(sum( D(p(BEG),p(i)) ));
      delta = delta - span * sum(sum( D(p(END),p(i)) ));

      edge  = (span+1) - 2*(1:span);
      delta = delta + edge*D(p(MID),p(i));

    else
      BEG = (1:j-1)';
      MID = (j:i-1)';
      END = (i+1:n)';

      delta = sum(sum( D(p(BEG),p(MID)) )) - sum(sum( D(p(END),p(MID)) ));

      delta = delta - span * sum(sum( D(p(BEG),p(i)) ));
      delta = delta + span * sum(sum( D(p(i),p(END)) ));

      edge  = 2*(1:span) - (span+1);
      delta = delta + edge*D(p(MID),p(i));

    end
    
  else
    n     = size(D,1);
    delta = 0;
    if j>i
      SPAN = j-i;
      L = (i+1:j)';
      M = (j+1:n)';
      delta = delta + sum(sum( D(p(M),p(L)) ));
      M = (1:i-1)';
      delta = delta - sum(sum( D(p(M),p(L)) ));
      M = (1:i-1)';
      delta = delta + SPAN * sum(sum( D(p(M),p(i)) ));
      M = (j+1:n)';
      delta = delta - SPAN * sum(sum( D(p(i),p(M)) ));

      SPAN2 = SPAN + 1;
      for m = i+1:j
        SPAN2 = SPAN2-2;
        delta = delta + SPAN2 * D(p(i),p(m));
      end
    else
      SPAN = i-j;
      L = (j:i-1)';
      M = (i+1:n)';
      delta = delta - sum(sum( D(p(M),p(L)) ));
      M = (1:j-1);
      delta = delta + sum(sum( D(p(M),p(L)) ));
      M = (1:j-1)';
      delta = delta - SPAN * sum(sum( D(p(M),p(i)) ));
      M = (i+1:n)';
      delta = delta + SPAN * sum(sum( D(p(i),p(M)) ));

      SPAN2 = SPAN + 1;
      for M = j:i:-1
        SPAN2 = SPAN2 - 2;
        delta = delta - SPAN2 * D(p(i),p(M));
      end
    end
  end

end

function [i,j] = random_nodes(min,max,order)
  %pick two random integers between (and including) min and max
  %sort the values if positive order
  pick = (min:max)';
  p    = randperm(max-min+1);
  p    = p(1:2);
  temp = pick(p);
  if order > 0
    temp = sort(temp);
  end
  i = temp(1);
  j = temp(2);
end

%--------------------------------------------------------------------------------------------------












