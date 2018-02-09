%AXIS algorithm (maybe known as Reciprocal Averaging Method)
%Independntly developed in three languages:
%Goldmann (1973) in Germany, Axis Wilkinson (1974) in U.K., and Legoux (1980) in France.
%GOLDMANN, K. 1973. "Zwei Methoden chronologisher Gruppierung", Acta Praehist et ArchaeoL, 3: 1-34.
%WILKINSON, E. M. 1974. "Techniques of Data Analysis-Seriation theory", Archaeo-physica, 5: 7-142.
%LEGOUX, R. 1980. in Perin, R, (ed.), La datation des tombes Mérovingiennes. Librairie Droz, Genève.
%
%Example (Fig.1 from Brower's paper):
%         X = [0 1 0 1 0; 1 1 0 1 0; 1 1 0 0 0; 0 0 1 0 1; 0 0 0 1 0; 1 1 0 0 0 ; 0 0 1 0 1; 0 0 1 1 0; 1 1 0 1 0];
%
%X       : an mxn matrix with the two-mode occurence or abundance data (should be non-negative, scale doesn't matter)
%'method': 'standard' or 'stabilised'
%perm    : the permutation vectors

%AJ Brockmeier & JY Goulermas, 2016

function [perm, details] = axis(X, varargin)

  if mod(length(varargin),2)
    error('dma:axis', 'missing parameters')
  end
  
  method    = 'standard'; %default values
  max_iters = 1000;

  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'method'
        method = varargin{i+1};
      case 'max_iters'
        max_iters = varargin{i+1};
      otherwise
        error('dma:axis', 'unknown parameter: %s', varargin{i})
    end
  end
  
  switch find([strncmpi(method, {'standard','improved'}, length(method)),1], 1)
    case 1
      [perm, details] = AxisStandard(X, max_iters);
    case 2
      [perm, details] = AxisImproved(X, max_iters);
    otherwise
      error('dma:axis', 'unknown method: %s', method)
  end
  
end
%--------------------------------------------------------------------------------------------------
%The stabilised AXIS method (by AJ Brockmeier), which often shows better 2-mode seriation scores
function [perm, details] = AxisImproved(X, max_iters)
  [m,n] = size(X);
  Q     = bsxfun(@times, X, 1./sum(X,1)); %Distribution of feature across instances
  P     = bsxfun(@times, X, 1./sum(X,2)); %Multinomial distribution of feature given instance
  x     = (1:m)';                         %location of instances

  for j = 1 : max_iters
    [~,p] = sort( P * (Q'*x) );
    ip    = dma.invperm(p);
    if isequal(x, ip)
      break;
    end
    x = ip;
  end
  
  perm.row      = dma.invperm(x);
  [~,perm.col]  = sort(Q'*x);
  details.steps = j;
end
%--------------------------------------------------------------------------------------------------
%The classic AXIS method, which can exhibit some instability
function [perm, details] = AxisStandard(X, max_iters)
  [m,n] = size(X);
  rows  = (1:m)';
  cols  = (1:n)';
  R     = 1./sum(X,2);
  C     = 1./sum(X,1)';
  stop  = false;
  iters = 0;
  while ~stop && iters < max_iters
    iters    = iters + 1;
    prev_sol = [rows;cols];
    %[~,rows] = sort(X*cols .* R); %buggy
    [~,rows] = sort(X *dma.invperm(cols) .* R); %or [~,rr] = sort( X(rows,:)*cols .* R(rows) );                rows = rows(rr);
    [~,cols] = sort(X'*dma.invperm(rows) .* C); %or [~,cc] = sort( X(:,cols)'*dma.invperm(rows) .* C(cols) );  cols = cols(cc);
    stop     = isequal(prev_sol, [rows;cols]);
  end
  
  perm.row      = rows;
  perm.col      = cols;
  details.iters = iters;
end
%--------------------------------------------------------------------------------------------------























