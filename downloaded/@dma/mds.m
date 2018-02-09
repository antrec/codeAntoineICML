%From paper: Hahsler et al., Getting things in order, JSS, 2008
%
%Orders the rows using the one component MDS; uses Matlab's cmdscale() & mdscale() functions.
%D        : an nxn distance matrix (see Matlab's MDS routines for further restrictions).
%'method' : stress criterion: 'classical', 'nonclassical', 'nonmetric'.
%'dims'   : dimensions to run non-classical MDS (dims=1 always gives better seriation quality).
%perm     : the permutation vector

%JY Goulermas, 2013

function [perm, details] = mds(D, varargin)

  if mod(length(varargin),2)
    error('dma:mds', 'missing parameters')
  end

  method = 'classical'; %default values
  dims   = 1;
  start  = 'cmdscale';
  
  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'method'
        method = varargin{i+1};
      case 'dims'
        dims = varargin{i+1};
      case 'start'
        start = varargin{i+1};
      otherwise
        error('dma:mds', ['unknown parameter: ', varargin{i}] )
    end
  end

  try
    switch lower(method)
      case 'classical'
        [Z, details.eigv] = cmdscale(D);
      case 'nonclassical'
        [Z, details.stress] = mdscale(D, dims, 'criterion', 'metricstress', 'start', start);
      case 'nonmetric'
          [Z, details.stress] = mdscale(D, dims, 'criterion', 'stress', 'start', start);
          %[Z, details.stress] = mdscale(D, dims, 'criterion', 'sstress', 'start', start, 'options', statset('display', 'iter') ); %experimental to avoid non-termination
    end%switch
    [~, perm] = sort(Z(:,1)); %size(Z,2) could be >1 for cmdscale(), and mdscale() if dims>1
  catch err
    if ( strcmp(err.identifier, 'stats:mdscale:ColocatedPoints') )   %nonclasiscal/nonmetric scaling fails when distances are zero
      [perm, details] = dma.mds(D, varargin{:}, 'start', 'random' ); %watch it here: even if 'start' is previously given as 'cmdscale', then 'random' is re-assigned later due to order!
    else
      rethrow(err);
    end
  end%try-catch
  
end