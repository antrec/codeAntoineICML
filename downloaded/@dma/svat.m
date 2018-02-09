%From paper: sVAT, scalable visual assessment of cluster tendency for large datasets, Pattern Recognition, Hathaway et al., 2006
%
%D         : a nxn dissimilarity matrix
%'clusters': an overestimate of true cluster number
%'sample'  : the new sample size (if <1, it is considered a fraction)
%perm      : the permutation vector (smaller than n)

%JY Goulermas, 2009

function [perm] = svat(D, varargin)

  clusters = 5;   %default values
  sample   = 0.6;
  
  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'clusters'
        clusters = varargin{i+1};
      case 'sample'
        sample = varargin{i+1};
      otherwise
        error('dma:svat', ['unknown parameter: ', varargin{i}] )
    end
  end

  n = size(D,1);
  if sample <= 1
    sample = sample * n;
  end
 
  %Step 1: select 'clusters' distinguished objects
  seed(1)   = 1;
  min_dists = D(1,:);
  for t = 2 : clusters
    min_dists  = min( min_dists, D(seed(t-1),:) ); %remember the closest points to the entire group of seeds
    %assert( isequal( min( D(seed,:), [], 1 ), min_dists ) )
    [~,seed(t)] = max(min_dists);
  end
  
  %Step 2: group objects with the nearest distinguished ones
  O = cell(1,clusters);
  for t = 1 : n
    [~, k] = min( D(seed,t) ); %for each datum, find closest seed
    O{k}   = [ O{k}, t ];
  end
  %assert( sum( cellfun(@length, O) ) == n )
  
  %Step 3: sample data
  U = [];
  for t = 1 : clusters
    new_size = ceil( sample * numel(O{t}) / n );
    idx      = randperm( numel(O{t}) );
    U        = [ U, O{t}( idx(1:new_size) ) ];
  end
  
  %Step 4: call VAT on the mxm small principal submatrix
  PS                  = D(U,U);
  [sample_perm, ~]    = dma.vat(PS);
  perm                = U(sample_perm); %map the VAT permutation of the sample to the original space
  assert( isequal( D(perm,perm), PS(sample_perm,sample_perm) ) )

end