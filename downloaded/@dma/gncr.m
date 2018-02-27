%From paper: Approximation Methods for Large Scale Object Sequencing.
%            X. Evangelopoulos, A.J. Brockmeier, T. Mu and J.Y. Goulermas,
%            (under submission)
%Input:
%   D      : an nxn symmetric dissimilarity matrix
%   method : 'twosum' to solve the 2-SUM and 'pshuber' to solve the Huberized 1-SUM
%   perm   : final permutation
%   details: struct with final 2-SUM/1-SUM scores, etc.
%
%J. Y. Goulermas & X. Evangelopoulos, 2018.


function [perm, details] = gncr(D, varargin)

  im = 1;
  while im <= length(varargin)-1 && ~strcmpi(varargin{im}, 'method')
    im = im + 1;
  end

  if im > length(varargin)-1
    error('dma:gncr', 'missing method specifier')
  else
    method            = varargin{im+1};
    varargin(im:im+1) = [];
  end

  switch lower(method)
    case 'twosum'
      [perm, details] = gncr_twosum(D, varargin{:});
    case 'pshuber'
      [perm, details] = gncr_pseudo_huber(D, varargin{:});        
    otherwise
        error('dma:gncr', 'unknown gncr type: %s', method )
  end

end
    