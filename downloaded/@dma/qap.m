%QAP-based seriation using various optimisers to solve the underlying QAP, given a seriation template.
%D   : an nxn distance matrix
%perm: the permutation vector

%JY Goulermas, 2013

function [perm, details] = qap(D, varargin)

  n         = size(D,1);
  method    = 'fw'; %default values
  max_iters = 999;
  W         = dma.template(n,'linear');
  
  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'method'
        method = varargin{i+1};
      case 'max_iters'
        max_iters = varargin{i+1};
      case 'w'
        W = varargin{i+1};
      otherwise
        error('dma:qap', 'unknown parameter: %s', varargin{i} )
    end
  end

  switch lower(method)
    case 'fw' %Frank-Wolfe approximation
      [perm, details] = dma.fw(D, W, max_iters);
      %Below is old caller using Vogelstein's sfw() code
      %[details.cost, perm, ~, details.iters,~, ~] = sfw(W, D, max_iters); %minimises sum(sum(W.*D(p,p)));
      
    case 'sa' %Simulated annealing
      [perm, details] = dma.sa(D, W, 'display_freq', 500, 'T_max', 500);

    case 'lopi' %Locally optimal pairwise interchange
      [perm, details] = dma.lopi(D, W, randperm(n));
      
    case 'fa' %Firefly algorithm, calling Yarpiz.com(c) code
      [perm, details] = yarpiz_fa(D, W, max_iters);
      
    case 'pso' %Particle Swarm Optimiser calling Yarpiz.com(c) code
      [perm, details] = yarpiz_pso(D, W, max_iters);

    case 'aco' %Ant Colony Optimiser calling Yarpiz.com(c) code
      [perm, details] = yarpiz_aco(D, W, max_iters);

    case 'pso2' %PSO calling Matlab
       options               = optimoptions(@particleswarm);
       options.Display       = 'iter';
       options.MaxIterations = max_iters;
       [s, ~, ~, details]    = particleswarm(@Objective, n, zeros(n,1), ones(n,1), options);
       [~,perm]              = sort(s);

    case 'ps' %Pattern Search calling Matlab's Global Optimisation toolbox
       options               = optimoptions(@patternsearch);
       options.Display       = 'iter';
       options.MaxIterations = max_iters;
       [s, ~, ~, details]    = patternsearch(@Objective, rand(n,1), [],[],[],[], zeros(n,1), ones(n,1), options); %needs the bounds!
       [~,perm]              = sort(s);

    case 'test2' 
       [perm, details] = dma.test2(D, W);
      
    otherwise
      error('dma:qap', 'unknown method: %s', method )
  end

  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  %Nested objective function where <P*D*P',W> is evaluated by
  %converting a continuous vector s to permutation matrix P. Used by some
  %methods made for continuous search.
  function [val] = Objective(s)
    [~,p] = sort(s);
    val   = sum(sum(D(p,p).*W));
  end
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
end





















