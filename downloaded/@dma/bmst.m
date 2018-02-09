%Bipartite MST, from Mu & Goulermas, ACAS paper, TPAMI, 2013
%
%X      : an mxn matrix with the two-mode data.
%perm   : the permutation vectors

%JY Goulermas, 2010

function [perm, details] = bmst(D)

  [m,n]     = size(D);
  in_rows  =  false(1,m) ;
  in_cols  =  false(1,n) ;
  perm_rows = [];
  perm_cols = [];  

  %initialise the row and column seeds as the most certain co-pair
  [~,idx]      = min( D(:) );
  [row,col]    = ind2sub( [m,n], idx );
  perm_rows(1) = row;
  perm_cols(1) = col;
  in_rows(row) = 1;
  in_cols(col) = 1;  

  while numel(perm_rows) + numel(perm_cols) < m + n
 
    val_c = min( min( D( in_rows,~in_cols) ) );
    val_r = min( min( D(~in_rows, in_cols) ) );
  
    if isempty(val_c), val_c = inf; end %one list may run out faster
    if isempty(val_r), val_r = inf; end
    
    if val_c < val_r %column node to enter the co-cluster
      [~,idx_c]      = find( D(in_rows,~in_cols) == val_c, 1, 'first' );
      indices        = find(~in_cols);
      idx_c          = indices(idx_c);
      perm_cols      = [perm_cols, idx_c];
      in_cols(idx_c) = 1;
    else %row node to enter
      [idx_r,idx_c]  = find( D(~in_rows,in_cols) == val_r, 1, 'first' );
      indices        = find(~in_rows);
      idx_r          = indices(idx_r);
      perm_rows      = [perm_rows, idx_r];
      in_rows(idx_r) = 1;
    end

  end
    
  perm.row  = perm_rows';
  perm.col  = perm_cols';
  details   = [];

end