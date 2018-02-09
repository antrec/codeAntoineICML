%From paper: Tsafrir et al., Sorting points into neighborhoods (SPIN): data analysis and visualisation by
%            ordering distance matrices, Bioinformatics, 2005
%
%D       : the nxn dissimilarity matrix
%'method': 'sts' or 'neighbourhood'
%'exact' : true iff the LAP is used for the Neighbourhood version of SPIN
%perm    : the resulting permutation vector

%JY Goulermas, 2011

function [perm, details] = spin(D, varargin)

  if mod(length(varargin),2)
    error('dma:spin', 'missing parameters')
  end
  
  method = 'sts'; %default values
  exact  = false;

  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'method'
        method = varargin{i+1};
      case 'exact'
        exact = varargin{i+1};
        if ~islogical(exact), error('dma:spin', 'logical parameter is required'), end
      otherwise
        error('dma:spin', 'unknown parameter: %s', varargin{i})
    end
  end
  
  switch find([strncmpi(method, {'sts','neighbourhood'}, length(method)),1], 1)
    case 1
      [perm, details.iters] = SideToSide(D);
    case 2
      [perm, details.iters] = Neighbourhood(D, exact);
    otherwise
      error('dma:spin', 'unknown method: %s', method )
  end
  
end

%--------------------------------------------------------------------------------------------------

%Side-to-side version
function [perm, iter] = SideToSide(D)
  n      = size(D,1);
  w      = (1:n)' - (n+1)/2; %implying W=w*w'
  P_prev = speye(n);
  stop   = false;
  iter   = 0;
  while ~stop
    iter   = iter + 1;
    s      = D * P_prev' * w;
    [~, p] = sort(s, 'descend');
    P      = dma.matperm(p);
    stop   = norm(P*s - P_prev*s, inf) < 1e-7;
    %imagesc(P*D*P'), obj = @(T) w'*T*D*T'*w; obj(P)-obj(P_prev), pause %%%debug
    P_prev = P;
  end
  perm = dma.matperm(P);
  %cost = w'*P*D*P'*w
end

%--------------------------------------------------------------------------------------------------

%Neighbourhood version
function [perm, iter] = Neighbourhood(D, exact_lap)
  n         = size(D,1);
  W0        = dma.template(n, 'squared');% -(r-c).^2;
  sigma     = 50; %smoothing factor
  resol     = max(D(:))/n; %to speedup the LAP (sometimes)
  P         = speye(n);
  cost_prev = inf;
  stop      = false;
  iter      = 0;
  while ~stop
    iter      = iter + 1;
    sigma     = max(0.2, sigma/2);   %reduce sigma gradually
    W         = exp(W0 / (n*sigma)); %I found that the smoothed maps from this gradual reduction is also achievable with a banded template and reducing the radius
    if ~exact_lap
      [~, ~, W] = SinkhornKnopp(W); %normalise to doubly stochastic (improves mainly the 'fast' case)
    end
    A         = D * P' * W;
    switch exact_lap
      case true                     %lapjv(A) minimizes trace(A*P') = trace(A(:,p)), where P=dma.matperm(p),
        [p,cost] = lapjv(A, resol); %or trace(P*A), where P'=dma.matperm(p)

      case false %approximate LAP
        [~,q]  = min(A, [], 2); %get minimum of each row and then break ties
        cur    = 1;
        labels = unique(q);
        p      = NaN(n,1);
        for k = labels'
          pos    = find(q==k);
          p(pos) = (cur : cur+numel(pos)-1)';
          cur    = cur + numel(pos);
        end%for
        cost = trace( A(:,p) );
    end%switch
    P    = dma.matperm(p)';
    stop = abs(cost-cost_prev) < 1e-7;
    if iter > 100
      stop = true; %iteration limit is needed for the inexact mode, as it may oscillate and never terminate!
      warning('dma:spin', 'iteration exceeded (%i total)', iter)
    end
    %iter, imagesc(P*D*P'); abs(cost-cost_prev), cost, pause %%%debug
    cost_prev = cost;
  end
  perm = dma.matperm(P);
end

%--------------------------------------------------------------------------------------------------
























