%Flips permutations that represents object sequences to make them more comparable.
%
%If 'type'='weighted', then:
%p can be a *column* vector, or an nxb matrix of b permutation vectors (the operation is
%applied column-wise). Each column is mirrored to make tail heavier.
%
%If 'type'='signed', then:
%p can only be nxb, with b>=2. In this case, the last b-1 columns are
%flipped where needed to make Spearman's rho of the same sign as the first column.

%JY Goulermas & A Kostopoulos, 2013

function [p] = flip(p, type)

  if ~exist('type', 'var')
      type = 'weighted';
  end

  [n,b] = size(p);
  
  switch lower(type)
    
    case 'weighted'  
      t = floor(n/2); %k could be lower, but larger will not change the result
      if isvector(p)
        if sum( p(1:t) ) > sum( p(end-t+1:end) ) %compare first t, with last t elements
          p = flipud(p);
        end
      else
        for k = 1 : b
          p(:,k) = dma.flip( p(:,k) );
        end
      end
      
    case 'signed'
      if b < 2,
        error('this type requires multiple vectors')
      end
      Inv = @(x) (x);
      ip1 = Inv( p(:,1) );
      for k = 2 : b
        sgn = dma.distperm( ip1, Inv( p(:,k) ), 'type', 'spearman_norm');
        if sgn < 0
          p(:,k) = flipud( p(:,k) );
        end
      end%for
      
    otherwise
      error('dma:flip', 'unknown permutation flipping type: %s', type )
  end%switch
  
end