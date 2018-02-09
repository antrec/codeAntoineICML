%From papers: OLO   : Bar-Joseph et al., Fast optimal leaf ordering, Bioinformatics, 2001
%             HC+R2E: Tien et al., Methods for simultaneously identifying coherent local clusters..., BMC Bioinformatics, 2008
%
%D         : an nxn distance matrix.
%'method'  : as in Matlab's linkage() for the fusion of nodes during hierarchical clustering: 'single', 'average', 'complete', etc.
%'ordering': how to order the terminal nodes: 'none' to retain original tree terminal order, 'olo' for optimal leaf ordering, 'external' to adopt external reference
%'ext_perm': the external reference vector when 'method' is 'external'
%perm      : the permutation vector

%JY Goulermas, 2013

function [perm, details] = hc(D, varargin)

  if mod(length(varargin),2)
    error('dma:hc', 'missing parameters')
  end

  method   = 'average'; %default values
  ordering = 'olo';
  
  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'method'
        method = varargin{i+1};
      case 'ordering'
        ordering = varargin{i+1};
      case 'ext_perm'
        ext_perm = varargin{i+1};
      otherwise
        error('dma:hc', ['unknown parameter: ', varargin{i}] )
    end
  end

  n    = size(D,1);
  tree = linkage(D, method);
  %clf, dendrogram(tree, 0, 'colorthreshold','default'); %if you like to see a tree version

  switch lower(ordering)
    case 'none'
      perm = GetLeaves(tree, n-1 ); %[~,~,perm] = dendrogram(tree,0); %another order with output
    case 'olo'
      perm = optimalleaforder(tree, D)';
    case 'external'
      perm = GetAdjustedLeaves(tree, n-1, dma.invperm(ext_perm) );
  end

  details.tree = tree;
  
end

%--------------------------------------------------------------------------------------------------

%Find leaves of a node with row index j (based on the way Matlab's linkage() encodes the tree)
function [leaves] = GetLeaves(T, j)

  n = size(T,1)+1;
  
  %linkage() encode a new node of index i at row n+i
  if T(j,1) <= n, left = T(j,1);
  else            left = GetLeaves(T, T(j,1)-n);
  end

  if T(j,2) <= n, right = T(j,2);
  else            right = GetLeaves(T, T(j,2)-n);
  end

  leaves = [left,right]; %any arbitrary order is fine here

end

%--------------------------------------------------------------------------------------------------

%Finds the leaves as the GetLeaves() function, but orders them according to an external reference permutation vector.
%P must be the inverse permutation, as the external reference wants i in the P(i)th position.
function [leaves] = GetAdjustedLeaves(T, j, P)

  n = size(T,1) + 1;
  
  if T(j,1) <= n, left = T(j,1);
  else            left = GetAdjustedLeaves(T, T(j,1)-n, P);
  end

  if T(j,2) <= n, right = T(j,2);
  else            right = GetAdjustedLeaves(T, T(j,2)-n, P);
  end

  ml = mean( P(left ) ); %Test here is the one from the HC_R2E paper.
  mr = mean( P(right));
  if ml <= mr, leaves = [left,  right];
  else         leaves = [right, left ];
  end

end

%--------------------------------------------------------------------------------------------------











































