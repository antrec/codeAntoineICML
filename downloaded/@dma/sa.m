%A Simulated Annealing implementation to minimise single- or multi-commodity QAPs expressed as: Sum_i trace( P* D_i *P' * F_i ).
%D, F : the distance and the flow matrices. These are either matrices for standard QAP problems, or cells of equal length
%       for multi-commodity ones.
%       Both D & F must be *symmetric* and (at least one of them) with *zero* diagonals (otherwise the deltas ar wrong!!!).
%       The operators are implemented to be (almost) as efficient as possible in terms of calculating the objective
%       difference for the changed part only.
%perm    : the optimising permutation
%details : history information
%p0      : the initial permutation
%T_max   : maximum temperature
%T_min   : freezing temperature to control termination
%T_length: iterations taken at each single temperature
%alpha   : annealing schedule parameter

%JY Goulermas, 2014

function [perm, details] = sa(D, F, varargin)
  if ~iscell(D)
    D = {D};
    F = {F};
  end
  L = numel(D);
  n = size(D{1},1);

  %default values
  p0           = randperm(n)'; %starting solution
  T_max        = 500;          %starting temperature (better do this adaptive)
  T_length     = 200;          %number of iterations taken at each temperature 
  T_min        = 1E-4;         %terminate when this is reached
  alpha        = 0.97;         %annealing parameter
  k            = 1.0;          %Boltzmann constant
  display_freq = 300;          %screen update frequency

  if mod(length(varargin),2)
    error('dma:sa', 'missing parameters')
  end

  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case {'p0','p', 'p_init'}
        p0 = varargin{i+1};
      case 't_max'
        T_max = varargin{i+1};
      case 't_min'
        T_min = varargin{i+1};
      case 't_length'
        T_length = varargin{i+1};
      case 'alpha'
        alpha = varargin{i+1};
      case 'display_freq'
        display_freq = varargin{i+1};
     otherwise
        error('dma:sa', 'unknown parameter: %s', varargin{i})
    end
  end

  if isequal(lower(T_max), 'auto')
    T_max = EstimateStartingTemperature(D,F);
  end
  
  time      = cputime;
  t         = 1;
  T         = T_max;
  p         = p0; %current solution and its cost
  p_best    = p0;
  cost      = Cost(D, F, p0); %initial cost calculated from scratch
  cost_best = cost;
  better    = 0; %counters
  worse     = 0;
  accept    = 0;
  stop      = false;
  while ~stop
    option = randi(7);
    switch option
      case {1,2,3} %Swap operator (this seems to be more effective than the others)
        [i,j]      = GetRandomIndices(1, n);
        idx        = [1:i-1, i+1:j-1, j+1:n]; %indices affected by the swap (simplified due to symmetry and zero diagonals)
        delta_swap = @(D,F) 2*( -D(p(i),p(idx)) +D(p(j),p(idx)) ) * ( F(i,idx) - F(j,idx) )'; %=new - old cost
        delta = 0;
        for k = 1 : L
          delta = delta + delta_swap(D{k},F{k});
        end
        p_new = p; p_new(i) = p(j); p_new(j) = p(i);
        
      case {4,5} %Inversion operator
        [i,j] = GetRandomIndices2(1, n, 5); %shorter inversion segment lengths give better results (perhaps adapt with t)
        idx1  = 1:i-1;
        idxN  = i:j;
        idxR  = j:-1:i;
        idx3  = j+1:n;
        delta_invert = @(D,F) -2*sum(sum( D(p(idxN),p(idx1)).*F(idxN,idx1) )) -sum(sum( D(p(idxN),p(idxN)).*F(idxN,idxN) )) -2*sum(sum( D(p(idxN),p(idx3)).*F(idxN,idx3) )) ...
                              +2*sum(sum( D(p(idxR),p(idx1)).*F(idxN,idx1) )) +sum(sum( D(p(idxR),p(idxR)).*F(idxN,idxN) )) +2*sum(sum( D(p(idxR),p(idx3)).*F(idxN,idx3) ));
        delta = 0;
        for k = 1 : L
          delta = delta + delta_invert(D{k},F{k});
        end
        p_new = p([idx1, idxR, idx3]);
        
      case {6,7} %Shift operator
        [i,j] = GetRandomIndices(1, n);
        idx1  = 1:i-1;
        idxN  = i:j;
        idx3  = j+1:n;
        if rand > 0.5 %left or right shift by 1 within [i,j]
          idxS  = [i+1:j, i];
        else
          idxS = [j, i:j-1];
        end
        delta_shift = @(D,F) -2*sum(sum( D(p(idxN),p(idx1)).*F(idxN,idx1) )) -sum(sum( D(p(idxN),p(idxN)).*F(idxN,idxN) )) -2*sum(sum( D(p(idxN),p(idx3)).*F(idxN,idx3) )) ...
                             +2*sum(sum( D(p(idxS),p(idx1)).*F(idxN,idx1) )) +sum(sum( D(p(idxS),p(idxS)).*F(idxN,idxN) )) +2*sum(sum( D(p(idxS),p(idx3)).*F(idxN,idx3) ));
        delta = 0;
        for k = 1 : L
          delta = delta + delta_shift(D{k},F{k});
        end
        p_new = p([idx1, idxS, idx3]);
    end%switch
    %%%assert( abs( cost + delta - Cost(D,F,p_new) ) < 1e-5, 'delta-swap, -invert or -shift is wrong!? (operator:%d)', option)

    %Jump state
    if delta <= 0 %better solution
      p    = p_new;
      cost = cost + delta;
      better = better + 1;
      if cost < cost_best
        cost_best = cost;
        p_best    = p_new;
      end
    else %worse solution
      rpE   = exp( - delta / (k*T) );
      worse = worse + 1;
      if rpE > rand
        p      = p_new;
        cost   = cost + delta;
        accept = accept + 1;
      end
    end
    
    if ~mod(t,display_freq) %screen update
      fprintf('\n Iteration: %d, at temperature: %2.3f, with minimising cost: %2.5f', t, T, cost);
    end
    
    if ~mod(t,T_length) %update temperature interval
      T = T * alpha;
    end
    details.cost_best(t) = cost_best;
    details.cost(t)      = cost;
    details.T(t)         = T;
    t                    = t + 1;
    stop                 = T < T_min;
  end%main SA loop
  
  perm          = p_best;
  details.time  = cputime - time;
  details.jumps = [better, worse, accept];
  %%figure, plot(1:t, details.cost, 'r', 1:t, details.cost_best, 'k')
end
%- - - - - - - - - - - - - - - - - - - - - - - - - - -
%Returns the overal multi-commodity QAP cost
function [cost] = Cost(D,F,p)
  cost = 0;
  for k = 1 : numel(D)
    cost = cost + sum(sum( D{k}(p,p) .* F{k} ));
  end
end
%- - - - - - - - - - - - - - - - - - - - - - - - - - -
%Get two random integers i & j, within [lb,ub] and with i < j.
function [i,j] = GetRandomIndices(lb, ub)
  i = randi([lb,ub-1]);
  j = randi([i+1,ub]);
end
%- - - -
%As above but also with abs(i-j) <= max_dst.
function [i,j] = GetRandomIndices2(lb, ub, max_dst)
  i = randi([lb,ub-1]);
  j = randi([i+1,min(ub,i+max_dst)]);
end
%- - - - - - - - - - - - - - - - - - - - - - - - - - -
%Heuristic for finding the maximum temperature
function [T_max] = EstimateStartingTemperature(D, F)
  n     = size(D{1},1);
  p     = randperm(n);
  T_max = -inf;
  for k = 1 : 1E4
    [i,j]      = GetRandomIndices(1, n);
    idx        = [1:i-1, i+1:j-1, j+1:n]; %indices affected by the swap (simplified due to symmetry and zero diagonals)
    delta_swap = @(D,F) 2*( -D(p(i),p(idx)) +D(p(j),p(idx)) ) * ( F(i,idx) - F(j,idx) )'; %=new - old cost
    delta = 0;
    for k = 1 : numel(D)
      delta = delta + delta_swap(D{k},F{k});
    end
    T_max = max(T_max, delta); %keep the worse
  end
end
%- - - - - - - - - - - - - - - - - - - - - - - - - - -

















