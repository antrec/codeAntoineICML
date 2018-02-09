%From paper: Alter et al., Singular value decomposition for genome-wide expression data processing and modeling, PNAS, 2000
%
%Uses SVD-based projections to order the angular arrangement of the concept-samples and features
%X   : an mxn matrix with the two-mode data, or a square one-mode matrix.
%'mode'  : 1 or 2 to specify the operating mode
%'subset': the indices of the two singular vectors to be used for projection
%'gap'   : false for standard algorithm, true to order from the largest gap
%perm    : the permutation vector(s)

%JY Goulermas, 2012

function [perm] = svd(X, varargin)

  if mod(length(varargin),2)
    error('dma:svd', 'missing parameters')
  end

  subset = [1,2]; %default values
  gap    = false;
  mode   = 2;
  
  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'subset'
        subset = varargin{i+1};
        if numel(subset) ~= 2
          error('dma:svd', 'two subset indices are required')
        end
      case 'mode'
        mode = varargin{i+1};
      case 'gap'
        gap = varargin{i+1};
      otherwise
        error('dma:svd', ['unknown parameter: ', varargin{i}] )
    end
  end

  [U,~,V] = svd(X,0);

  Y = X * V(:,subset);
  a = atan2( Y(:,2), Y(:,1) );
  if ~gap, [~,perm] = sort(a);
  else     perm     = CircularSort(a);
  end

  if mode == 2
    Y = X' * U(:,subset);
    a = atan2( Y(:,2), Y(:,1) );
    if ~gap, [~,perm2] = sort(a);
    else     perm2     = CircularSort(a);
    end

    perm = struct('row', perm, ...
                  'col', perm2 ...
                 );
  end%mode

end
