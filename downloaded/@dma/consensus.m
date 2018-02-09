%Estimates a new ranking, which is 'as close as possible' to b given rankings, using various techniques.
%
%P       : the nxb matrix of b ranking/permutation vectors in columns
%w       : a b-length vector of non-negative weights with sum 1, to weight the rankings. If w==[] then equal weights are assumed.
%'method': the type of method to be invoked, followed by the name of the technique, such as: 'borda', 'cmedian', etc.

%A Kostopoulos & JY Goulermas, 2013

function [perm, details] = consensus(P, w, varargin)

  if mod(length(varargin),2)
    error('dma:consensus', 'missing parameters')
  end

  [n,b]           = size(P);
  method          = 'borda'; %default values
  ppm_algo        = 'spectral';
  ppm_ga_args     = {'generations', 500, 'pop_size', 200, 'tol_fun', 1e-7};
  kendall_algo    = 'ga';
  kendall_ga_args = {'generations', 500, 'pop_size', 200, 'tol_fun', 1e-7};
  
  if isequal(w,[]), w = ones(b,1)/b;
  else              w = w(:) / sum(w);
  end

  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'method'
        method = varargin{i+1};
      case 'ppm_algo'
        ppm_algo = varargin{i+1};
      case 'ppm_ga_args'
        ppm_ga_args = varargin{i+1};
      case 'kendall_algo'
        kendall_algo = varargin{i+1};
      case 'kendall_ga_args'
        kendall_ga_args = varargin{i+1};

      otherwise
        error('dma:consensus', 'unknown parameter: %s', varargin{i} )
    end
  end

  switch lower(method)

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %From paper: Borda (1781) Memoire sur les electrions au scrutin. Histoire de l'Academie Royale des Sciences
    case 'borda'
      [~, I]    = sort(P);
      scores    = (n+1-I) * w; %Borda scores (number of wins of each object)
      [~, perm] = sort(scores, 'descend');

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %From paper: Condorcet (1785) Essai sur l'application de l'analyse a la probabilite des decisions rendues a la pluralite des voix.
    case 'condorcet'
      count = zeros(n); 
      for k = 1 : b %raw Condorcet counts
        count = count + w(k) * dma.relation( P(:,k) );
      end
      net_count = count - count'; %net Condorcet counts
      [~, perm] = sort( sum(sign(net_count),1)' ); 

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %From paper: Copeland (1951) A reasonable social welfare function. University of Michigan
    case 'copeland'
      count = zeros(n); 
      for k = 1 : b
        count = count + w(k) * dma.relation( P(:,k) );
      end
      net_count = count - count';
      net_count = ( sign(net_count) + 1 ) / 2;
      [~, perm] = sort( sum(net_count,1)' );

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Exact optimisation using Kendall or KS distance
    case 'cmedian'
      perm = CMedian(P, w);


    %-----------------------------------------------------------
    %new stuff below
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Minimises average Hamming distance, using LAP
    case 'hamming'
      C = zeros(n);
      for k = 1 : b
        C = C + w(k) * dma.matperm( P(:,k) ); %weighted sum of permutation matrices
      end
 
      [p, cost]    = lapjv(-C); %minimise trace(-C'*P) == trace(-C(:,p)) == cost, where P=dma.matperm(p)
      perm         = p';
      details.cost = cost;

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Minimises average Spearman distance, using LAP
    case 'spearman'
      Pinv = dma.invperm(P);         %inverse permutations are needed here, to use the position differences and not the object labels
      s    = sum(Pinv * diag(w), 2); %note: S*e=s
      
      if false %optimise underlying LAP (unecessary)
        E            = s * (1:n);
        [p, cost]    = lapjv(-E); %minimise trace(-E'*P) == trace(-E(:,p)) == cost, where P=dma.matperm(p)
        perm         = dma.invperm(p)';
        details.cost = cost;
      else
        [~,perm] = sort(s, 'ascend'); %maximise q'*s, so q needs to be ordered invperm(perm)
      end

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Minimises average Kendalldistance, using QAP
    case 'kendall'
      switch lower(kendall_algo)
        case 'mip'
          perm = CMedian(P, w); %cmedian here
        case 'ga'
          hybs{1}         = dma.invperm( dma.consensus(P, w, 'method', 'spearman')' ); %seed with this solution to be at least as good as that
          kendall_ga_args = {kendall_ga_args{:}, 'hybrids', hybs};
          o               = ones(n,1);
          B               = @(x) x*o' - o*x';
          C               = @(x) sign( B(x) );
          e               = (1:n)';
          D               = C(e);          
          N               = zeros(n);
          for k = 1 : b
            p = dma.invperm(P(:,k));
            N = N + w(k) * C(p);
          end
          o               = dma(D, 'ga');
          o.ga_args       = {'objective', 'qap', 'W', -N, kendall_ga_args{:} };
          perm            = dma.invperm( o.run().perm );
        otherwise
          error('dma:consensus', 'unknown kendall algorithm: %s', kendall_algo )
      end%switch

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Maximises the positional proximity between (i,j) pairs that are mostly preferred by P
    case 'ppm'
      K = dma.pdm(P, w);
      switch lower(ppm_algo)
        case 'spectral'
          o               = dma(K, 'spectral'); %K=1./(n-W) or K=exp(W-n) seems to be needed instead, for 'ding'
          o.spectral_args = {'method', 'barnard'};
          perm            = o.run().perm;
        case 'ga'
          hybs{1}     = dma.invperm( dma.consensus(P, w, 'method', 'ppm', 'ppm_algo', 'spectral')' ); %seed with the spectral solution
          ppm_ga_args = {ppm_ga_args{:}, 'hybrids', hybs};                                            %to be at least as good as that
          e           = (1:n)';
          D           = bsxfun(@minus, e, e').^2;
          o           = dma(D, 'ga');
          o.ga_args   = {'objective', 'qap', 'W', -K, ppm_ga_args{:} };
          perm        = dma.invperm( o.run().perm );
        otherwise
          error('dma:consensus', 'unknown ppm algorithm: %s', ppm_algo )
      end%switch
      details.K = K;

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    otherwise
      error('dma:consensus', 'unknown consensus method: %s', method )
  end

end

%--------------------------------------------------------------------------------------------------

%From: Grötschel, et al., A Cutting Plane Algorithm for the Linear Ordering Problem, Operations Research, 1984.
%      Hornik & Meyer, Deriving Consensus Rankings from Benchmarking Experiments, 2006.
function [p] = CMedian(P,w)

  [n,b] = size(P);
  nv    = n*(n-1)/2;         %number of variables
  nc    = n*(n-1)*(n-2) / 6; %number of constraints (half of)

  %Form objective function
  C = zeros(n);
  for i = 1 : b  
    Ri = dma.relation( P(:,i) );
    C  = C + 2*w(i)*Ri;
  end
  C = C - C';                %we want to minimise sum(sum( triu( (C-C').*R ) ))  
  c = C( triu(true(n),+1) ); %upper triangular of C in vector form

  %Form inequality constraints as Ax <= b
  % + pij + pjk - pik <= 1 for 1 <= i < j < k <= n
  % - pij - pjk + pik <= 0 for 1 <= i < j < k <= n
  T   = nchoosek(1:n,3);
  pij = f( T(:,1), T(:,2) );
  pjk = f( T(:,2), T(:,3) );
  pik = f( T(:,1), T(:,3) );
  A   = zeros(nc, nv);
  A(sub2ind( [nc nv], [1:nc]', pij )) = +1;
  A(sub2ind( [nc nv], [1:nc]', pjk )) = +1;
  A(sub2ind( [nc nv], [1:nc]', pik )) = -1;
  A = [ A; -A ];
  b = [ ones(nc,1); zeros(nc,1) ];
 
  %Solve: max(c'*x) s.t. Ax<=b
  options = optimoptions(@intlinprog,...
                         'LPMaxIter',  1E6*nv, ...
                         'Display',    'off');
  opt = intlinprog(-c, 1:nv, A, b, [], [], zeros(nv,1), ones(nv,1), options );

  %Create relation matrix R from vector optimal vector opt
  R = zeros(n);
  R( triu( true(n), +1 ) ) = opt;
  R = R + tril( 1-R', -1 );

  [~,p] = sort( sum(R,1)' );


  %Convert matrix upper triangular indices into vector form
  function x = f(i,j)
    x = (j-1).*(j-2)/2 + i;
  end
  
end

%--------------------------------------------------------------------------------------------------





































