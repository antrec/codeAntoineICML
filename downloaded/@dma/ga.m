%Use a genetic algorithm to optimise various seration scores and
%other objective functions based on templates.
%
%D    : an nxn distance matrix.
%perm : the permutation vector

%JY Goulermas, 2013

function [perm, details] = ga(D, varargin)

  if mod(length(varargin),2)
    error('dma:ga', 'missing parameters')
  end

  n           = size(D,1);
  objective   = 'tsp'; %default values
  generations = 800;
  pop_size    = 100;
  stall_gen   = 50;
  tol_fun     = 1e-6;
  W           = dma.template(n,'linear');
  hybrids     = 'off';

  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'objective'
        objective = varargin{i+1};
     case 'w'
        W = varargin{i+1};
     case 'generations'
        generations = varargin{i+1};
     case 'pop_size'
        pop_size = varargin{i+1};
     case 'tol_fun'
        tol_fun = varargin{i+1};
     case 'stall_gen'
        stall_gen = varargin{i+1};
     case 'hybrids'
        hybrids = varargin{i+1};
      otherwise
        error('dma:ga', 'unknown parameter: %s', varargin{i} )
    end
  end
  
  objective = @(x) ObjectiveFunction(x, objective, D, W);
  creation  = @(x,y,z) CreatePopulation(x, y, z, hybrids, D);
  options   = gaoptimset('PopulationType',    'custom', ...
                         'PopulationSize',    pop_size, ...
                         'PopInitRange',      [1;n], ...
                         'CreationFcn',       creation, ...
                         'InitialPopulation', [], ...
                         'FitnessScalingFcn', @fitscalingrank, ...
                         'EliteCount',        2, ...
                         'SelectionFcn',      @selectionstochunif, ... NB:selectionstochunif:best
                         'CrossoverFraction', 0.5, ...
                         'CrossoverFcn',      @Crossover, ...
                         'MutationFcn',       @Mutation, ...
                         'Generations',       generations, ...
                         'TolFun',            tol_fun', ...
                         'StallGenLimit',     stall_gen, ...
                         'Vectorized',        'off', ...
                         'PlotFcns',          [], ...
                         'display',           'iter'); %iter or off
  %options  = gaoptimset(options, 'PlotFcns', {@gaplotbestf,@PlotDistance});

  [opt, details.fval, details.reason, details.output] = ga(objective, n, options);
  perm = opt{:}';

end
  
%--------------------------------------------------------------------------------------------------

%The minimising objective for the GA, supporting different objective functions
function [obj] = ObjectiveFunction(p, objective, D, W)
  
  p = p{:};
  n = length(p);
  switch objective
    
    case 'tsp' %solves the standard TSP
      obj = 0;
      for i = 1 : n-1
        obj = obj + D( p(i), p(i+1) );
      end
      
    case 'qap' %solves the standard QAP using template W
      obj = sum(sum( D(p,p) .* W ));
      
    case 'arc' %anti-robinson compatibility
      arc = 0;
      D   = triu(D(p,p),1);
      for i = 2 : n-1
        diff_hor = triu( D(1:n-i,2:n+1-i) - D(1:n-i,i+1:n) );
        diff_ver = triu( D(i:n-1,i+1:n)   - D(1:n-i,i+1:n) );
        arc      = arc - sum(sum( sign(diff_hor) + sign(diff_ver) ));
      end
      obj = -arc;

    case 'are' %anti-robinson events violation
      are = 0;
      D   = triu(D(p,p),1);
      for i = 2 : n-1
        diff_hor = triu( D(1:n-i,2:n+1-i) - D(1:n-i,i+1:n) );
        diff_ver = triu( D(i:n-1,i+1:n)   - D(1:n-i,i+1:n) );
        are      = are + sum(sum( (diff_hor>0) + (diff_ver>0) ));
      end
      obj = +are;

    otherwise
      error('dma:ga:ObjectiveFunction', 'unknown GA objective formulation: %s', objective )

  end
end

%--------------------------------------------------------------------------------------------------

%Create the new population with random and hybrid individuals
function [pop] = CreatePopulation(n, ~, options, hybrids, D)

  pop_size = sum(options.PopulationSize);
  pop      = cell(pop_size,1);
  
  for i = 01 : pop_size %fully random population
    pop{i} = randperm(n); 
  end

  if isa(hybrids, 'char')
    if sum( strcmpi(hybrids, {'off', ''}) )
      %do nothing; keep random population
    elseif strcmpi(hybrids, 'seriation') %these seed solutions are good for seriating D
      pop{1}    = dma(D,'vat').run().perm';
      pop{2}    = dma(D,'r2e').run().perm';
      pop{3}    = dma(D,'pca').run().perm';
      pop{4}    = dma(D,'mds').run().perm';
      o         = dma(D, 'hc');
      o.hc_args = {'method', 'average', 'ordering', 'olo'};
      pop{5}    = o.run().perm';
    else
      error('dma:ga:CreatePopulation', 'unknown hybridisation string' )
    end
  elseif isa(hybrids, 'cell')
    for i = 1 : length(hybrids)
      pop{i} = hybrids{i};
    end
  else
    error('dma:ga:CreatePopulation', 'unknown hybridisation data structure' )
  end

end

%--------------------------------------------------------------------------------------------------

%Implement multiple mutation schemes for permutation chromosomes
function [children] = Mutation(parents, ~, n, ~, ~, ~, pop)

  children = pop(parents); %prepare copies
  
  for i = 1 : size(children,1)
    switch randi(5) %choose a random mutation scheme (biases experimentally set)
      case {1} %MULTIPLE GENE SWAPS
        k     = randi(max(1,round(n/30))); %number of gene swaps in total per child
        genes = [];
        while length(unique(genes)) ~= k*2 %to avoid redundant mutations and invalid children
          genes = randi(n, [k,2]);
        end
        children{i}(genes(:,1)) = pop{parents(i)}(genes(:,2));
        children{i}(genes(:,2)) = pop{parents(i)}(genes(:,1));

      case {2,3} %SEGMENT FLIP
        start                   = randi(n-1);
        stop                    = randi(n-start) + start;
        children{i}(start:stop) = pop{parents(i)}(stop:-1:start);

      case {4,5} %FORWARD/BACKWARD SEGMENT SHIFT
        limit  = round(n/5);                    %adaptive limit as mutations can be disruptive
        start  = randi(n-2)+1;
        stop   = min(n-1, start + randi(limit) ); %right segment edge limited
        if rand > 0.5                             %+ve shift
          sh_max = min(+limit, n-stop);           %longest shift also limited
          shift  = randi(sh_max);
          idx    = [1:start-1, stop+1:stop+shift, start:stop, stop+shift+1:n]; %trailing-covered-slid-remaining
        else
          sh_min = max(-limit, 1-start);
          shift  = -randi(-sh_min);
          idx    = [1:start+shift-1, start:stop, start+shift:start-1, stop+1:n]; %trailing-slid-covered-remaining
        end
        children{i} = pop{parents(i)}(idx);
    end%switch
    %assert(isequal(sort(children{i}),1:n))
  end%for

end

%--------------------------------------------------------------------------------------------------

%Call multiple crossover schemes for permutation chromosomes
function [children] = Crossover(parents, ~, ~, ~, ~, pop)

  m        = length(parents) / 2; %1 child only per parent pair
  children = cell(m,1);

  for i = 1 : m
    switch randi(3) %choose a random crossover scheme (biases experimentally set)
      case {1}
        children{i} = CrossoverPMX   ( pop{parents(2*i-1)}, pop{parents(2*i)} );
      case {2,3}
        children{i} = CrossoverOrder1( pop{parents(2*i-1)}, pop{parents(2*i)} );
    end%switch
  end%for

end

%--------------------------------------------------------------------------------------------------

%Creates a child p from permuation parents p1 & p2, using PMX crossover
function [p] = CrossoverPMX(p1, p2)

  n             = length(p1);
  p             = NaN(1,n);
  start         = randi(n-1);
  stop          = randi(n-start) + start;
  p(start:stop) = p1(start:stop);
  not_copied    = setdiff( p2(start:stop), p(start:stop) );

  for allele = not_copied
    index    = TracedIndex(allele);
    p(index) = allele;
  end

  safe    = isnan(p); %copy the remaining
  p(safe) = p2(safe);
  
  function [index] = TracedIndex(allele)
    value = p1(p2 == allele);
    index = find(p2 == value);
    if index >= start && index <= stop
      index = TracedIndex(value);
    end
  end

  %assert(isequal(sort(p),1:n))
end

%--------------------------------------------------------------------------------------------------

%Creates a child p from permutation parents p1 & p2, using Order1 crossover
function [p] = CrossoverOrder1(p1, p2)

  n             = length(p1);
  p             = NaN(1,n);
  start         = randi(n-1);
  stop          = randi(n-start) + start;
  p(start:stop) = p1(start:stop);

  [~,not_idx]             = setdiff(p2, p(start:stop));%setdiff row or col depending on Matlab version!
  not_idx                 = sort(not_idx);
  k                       = find([not_idx(:);inf] > stop, 1, 'first');
  p([stop+1:n,1:start-1]) = p2(not_idx([k:end,1:k-1]));
  
  %assert(isequal(sort(p),1:n))
end

%--------------------------------------------------------------------------------------------------

%To monitor genetic diversity (better permutation distance function, such as KS-distance is needed)
function [state] = PlotDistance(options, state, flag)

  max_pairs = 1000; %too slow to compare all m(m-1)/2 pairs; sampling needed
  m         = sum(options.PopulationSize);
  pairs     = nchoosek(1:m,2);
  p         = randperm(length(pairs));
  pairs     = pairs( p(1:min(max_pairs,length(pairs))), :);
  pop       = state.Population;
  dist      = 0;
  for p = 1 : length(pairs)
    dist = dist + abs( corr(pop{pairs(p,1)}', pop{pairs(p,2)}', 'type', 'kendall') );
  end
  dist = dist / length(pairs);

  switch flag %code below follows Mathworks's exact template
    case 'init'
      plotDist = plot(state.Generation, dist,'.');
      set(gca,'xlimmode','manual','zlimmode','manual', 'alimmode','manual')
      set(gca,'xlim',[1,options.Generations]);
      set(plotDist,'Tag','gaplotdistance');
      xlabel('Generation','interp','none');
      ylabel('Average Distance');
      title('Average Distance Between Individuals','interp','none')
    case 'iter'
      plotDist = findobj(get(gca,'Children'),'Tag','gaplotdistance');
      newX = [get(plotDist,'Xdata') state.Generation];
      newY = [get(plotDist,'Ydata') dist];
      set(plotDist,'Xdata',newX,'Ydata',newY);
  end

end

%--------------------------------------------------------------------------------------------------







