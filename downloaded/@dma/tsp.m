%TSP based seriation using various solvers. Formulation as in paper: Hahsler et al., Getting things in order, JSS, 2008
%
%D          : an nxn distance matrix.
%'starts'   : number of generated tours by the heuristic solver
%'precision': rounding precision for the concorde solver
%perm       : the permutation vector

%JY Goulermas, 2013

function [perm, details] = tsp(D, varargin)

  if mod(length(varargin),2)
    error('dma:tsp', 'missing parameters')
  end

  n         = size(D,1);
  method    = 'heuristic'; %default values
  starts    = n;
  precision = 2;

  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'method'
        method = varargin{i+1};
      case 'starts'
        starts = min(n, varargin{i+1});
      case 'precision'
        precision= varargin{i+1};
      otherwise
        error('dma:tsp', 'unknown parameter: %s', varargin{i} )
    end
  end

  D2 = [zeros(n+1,1), [zeros(1,n); D] ]; %insert a dummy node with zero connection costs to all others, to convert the Hamiltonian cycle to a path problem
  
  switch lower(method)
    case 'heuristic'
      [perm, details.cost] = tspsearch(D2, starts+1);

    case 'concorde' %external library (not included in this class!); compile & install yourself
      perm    = concorde('-x','-s','7481', ceil(D2*10^precision)) + 1;
      details = [];

    otherwise
      error('dma:tsp', 'unknown method: %s', method )
  end

  cut  = find(perm == 1);
  perm = perm([cut+1:end,1:cut-1]) -1;
  
end

%--------------------------------------------------------------------------------------------------

