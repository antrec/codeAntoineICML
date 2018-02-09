%From paper: Bezdek & Hathaway, VAT, a tool for visual assessment of (cluster) tendency, IJCNN, 2002
%
%D       : a nxn dissimilarity matrix
%dist_tol: if not zero, minimal distance value is relaxed to account for round-off errors
%perm    : the permutation vector

%JY Goulermas, 2010

function [perm, ties] = vat(D, dist_tol)

  if nargin == 1
    dist_tol = 0.0;
  end

  n    = size(D,1);
  perm = nan(n,1);
  ties = 0;
  
  %initialise
  [~,idx]   = max( D(:) ); %at least two ties here
  %[~,idx]  = min( D(:) ); %alternative for specific cases!
  [row,~]   = ind2sub([n,n], idx);
  perm(1)   = row;
  in        =  false(1,n);
  in(row)   = 1; %set of selected vertices
  out       = ~in;
  
  for t = 2 : n %assert( isequal(sort(perm(isfinite(perm))), find(in)') )
    min_val = min( min( D(in,out) ) );
    [near_in,near_out] = find( D(in,out) <= min_val + dist_tol ); %near_out = nearest vertex to add next
    if numel(near_out) > 1
      ties     = ties + 1;
      near_out = near_out( TieBreak('recent') );
    end
    idx           = find(out);     %map col to
    near_out      = idx(near_out); %original coordinates in D
    perm(t)       = near_out;
    in(near_out)  = 1;
    out(near_out) = 0;
  end
  
  %----------------------------------------------
  %TieBreak: 
  function pos = TieBreak(type)
    
    switch type
      case 'simple'
        pos = 1;
      case 'recent'
        in_idx   = find(in); %= sort(perm(1:t-1))
        near_in  = in_idx(near_in); %map to original coordinates
        [~, pos] = ismember( perm(1:t-1), near_in ); %find which in 'in' link to the nearest in 'out'
        pos      = pos( find(pos, 1, 'last') );
    end

  end
  %----------------------------------------------
  
end