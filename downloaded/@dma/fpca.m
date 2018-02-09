%Orders the rows using the first component from the PCA of D (for one-mode),
%and additionally the columns of D (for two-mode) with the first component of D'.
%D   : an nxm matrix with the two-mode data, or a square one-mode matrix.
%mode: 1 or 2 to specify the operating mode
%perm: the permutation vector(s)

%JY Goulermas, 2013

function [perm] = fpca(D, varargin)

  if mod(length(varargin),2)
    error('dma:fpca', 'missing parameters')
  end

  dims = 5; %default values
  mode = 1;

  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'dims'
        dims = varargin{i+1};
      case 'mode'
        mode = varargin{i+1};
      otherwise
        error('dma:fpca', ['unknown parameter: ', varargin{i}] )
    end
  end%for
  
  X       = bsxfun(@minus, D, mean(D,1) );
  [~,S,V] = svd(X,0);
  s       = diag(S);
  if dims == inf
    dims = sum( s > max(size(X)) * eps(max(s)) ); %use up to rank(X)
  end
  
  Y         = X * V(:,1:dims);
  [~, P]    = sort(Y,1); %permutations
  w         = s(1:dims); %weights  
  [perm, ~] = dma.consensus(P, w, 'method', 'ppm', 'ppm_algo', 'spectral');

  
  if mode == 2
    perm = struct('row', perm, ...
                  'col', dma.fpca(X', 'dims', dims, 'mode', 1) ...
                 );
  end

end
