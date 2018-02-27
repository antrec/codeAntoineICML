%Distance Matrix Analysis (DMA) class.
%Includes various algorithms, utilities and measures for one-mode and two-mode seriation,
%and processing of permutation information
%
%This is version: ver.2.9.4 (10-02-2018)

%JY Goulermas, 2010

classdef dma
  
  properties
    D             = [];    %nxm (one- or two-mode) dissimilarity matrix
    algo          = 'vat'; %main method to run by default
    mode          = 1;     %for 1-mode or 2-mode analysis; used by some methods only, when there is option for both cases, but ignored for others
    revat_sample  = 0.6;   %new sample size (if <1, is taken as a fraction)
    svat_args     = {};    %'clusters' for cluster over-estimate, 'sample' for new sample size (<1 means as a fraction)
    spin_args     = {};    %'method' is: 'sts' or 'neighbourhood', 'exact' is Boolean
    axis_args     = {};    %'method' is: 'standard' or 'improved', 'max_iters' maximum number of iterations
    hc_args       = {};    %'method' needs a method value as in 'linkage()', 'ordering' a value from {'none', 'olo', 'external'}, 'ext_perm' an external permutation vector
    mds_args      = {};    %'method' is: 'classical', 'nonclassical' or 'nonmetric'
    svd_args      = {};    %'subset' to set the two indices of sing. vectors, 'mode' and  'gap' true iff we order from largest gap
    fpca_args     = {};    %'dims' to set the number of components to use (inf for all up to the rank)
    tsp_args      = {};    %'starts' number of path start attempts in the heuristic solver, 'method' for 'heuristic' or 'concorde'
    qap_args      = {};    %'method' for 'sa' 'lopi','fw','fa','pso',pso2', or 'ps'. 'max_iters' for iterations. 'W' for a minimising seriation template.
    spectral_args = {};    %'method' is the name of the method, 'ext_perm' an external permutation vector, 'ext_mix' its mixture coefficient
    gncr_args     = {};    %'method' can be 'twosum' or 'pshuber', 'mu_unit' the mu parameter, 'gamma' the continuation rate, 'delta' for pshuber
    ga_args       = {};    %'objective' is one of: 'tsp', 'qap', era, arc', etc., W is the template for the 'qap' case; see corresponding function for more
    arsa_args     = {};    %simulated annealing method (see function for parameters)
    test1_args    = {};    %for testing - do NOT use!
  end
  
  properties (SetAccess = private, GetAccess = public)
    n                = 0;   %rows of D
    m                = 0;   %columns in D
    perm             = [];  %resulting permutation vector
    Q                = [];  %permuted version of D
    scores           = [];  %seriation measures of D & Q
    vat_ties         = NaN;
    revat_details    = [];
    covat_details    = [];
    r2e_details      = [];
    hc_details       = [];
    spin_details     = [];
    mds_details      = [];
    tsp_details      = [];
    qap_details      = [];
    spectral_details = [];
    gncr_details     = [];
    ga_details       = [];
    axis_details     = [];
    test1_details    = [];
    
  end
  
  properties (Constant = true)
    algos    = {'vat', 'revat', 'svat', 'covat', 'covat2', ...
                'corder', 'pca', 'svd', 'mds', 'bea', 'r2e', 'spin', ...
                'hc', 'tsp', 'spectral', 'gncr', 'qap', 'ga', 'arsa', 'fpca', ...
                'bmst', 'test1', 'axis'
               };
    dist_tol = 0.0; %used by vat to set tolerance
  end
    
  methods
    
    %ctor
    function obj = dma(x, algo)
      if ~nargin, return, end
      if isa(x,'dma'), obj.D = x.D;
      else             obj.D = x;
      end
      
      if exist('algo', 'var')
        obj.algo = algo;  
      end
      
    end

    %setter
    function obj = set.D(obj,D)
      obj.D         = D;
      [obj.n,obj.m] = size(D);
      obj.perm      = [];
      obj.Q         = [];
      obj.scores    = [];
      
      if obj.n ~= obj.m, obj.mode = 2;
      else               obj.mode = 1;
      end
    end

    %setter
    function obj = set.mode(obj,m)
      if ~ismember(m, [1,2])
        error('dma:set:mode', 'mode of %d is not supported', m)
      else
        obj.mode = m;
      end
      if obj.m ~= obj.n && obj.mode == 1
          error('dma:set:mode', '1-mode analysis for rectangular matrices not allowed')
      end
    end
    
    %setter
    function obj = set.algo(obj,algo)
      pos = strcmpi(algo, obj.algos);
      if sum(pos)
        obj.algo = algo;
      else
        error('dma:set:algo', 'unknown algorithm: %s', algo)
      end
      
    end
    
    %getter
    function value = get.scores(obj)
      value = obj.scores;
    end

    %analyse
    function obj = analyse(obj)
      obj.scores.D = dma.measures(obj.D);
      obj.scores.Q = dma.measures(obj.Q);
    end
    
    %subsref
    function val = subsref(obj, s)
      switch s(1).type
        case '.',  switch s(1).subs
                     case 'dims'
                       val = [obj.n, obj.m]';
                     case {'p', 'P'}
                       val = obj.perm;
                     case 'diff'
                       struct2mat = @(x) cell2mat(struct2cell(x));
                       val        = struct2mat(obj.scores.D) - struct2mat(obj.scores.Q);
                       otherwise
                         val = builtin('subsref', obj, s);
                   end
        otherwise, error('dma:subsref', 'unimplemented access route')
      end
    end

    %display
    function obj = graph(obj)
       fprintf('Input matrix D (%gx%g)\n\n', obj.n, obj.m);
       clf
       if obj.m == obj.n, axis_sc = 'square'; else axis_sc='equal'; end
       subplot(121), imagesc(obj.D), axis(axis_sc), title('D')
       subplot(122), imagesc(obj.Q), axis(axis_sc), title('Q')
    end
        
    %display for revat
    function obj = plot_revat_profiles(obj)
      figure
      piv  = obj.revat_details.pivot;
      prof = obj.revat_details.profile;
      for i = 1 : numel(piv)
       subplot( 2, ceil(numel(piv)/2), i )
       bar( prof(i,:) )
       axis tight
      end
    end
    
  end%methods
  
  %statics
  methods (Access = public, Static = true)
    %algorithms:
    [x,y] = vat      (x,y)
    [x,y] = revat    (x,y)
    [x]   = svat     (x,varargin)
    [x,y] = dbe      (x,y,z)
    [x,y] = covat    (x)
    [x,y] = covat2   (x)
    [x,y] = axis     (x,varargin)
    [x,y] = bmst     (x)
    [x,y] = r2e      (x,y)
    [x,y] = spin     (x,varargin)
    [x]   = corder   (x)
    [x]   = pca      (x,y)
    [x]   = svd      (x,varargin)
    [x]   = fpca     (x,varargin)
    [x,y] = mds      (x,varargin)
    [x]   = bea      (x,y)
    [x]   = bea_naive(x,y) %to be removed
    [x,y] = hc       (x,varargin)
    [x,y] = tsp      (x,varargin)
    [x,y] = qap      (x,varargin)
    [x,y] = spectral (x,varargin)
    [x,y] = gncr     (x,varargin)
    [x,y] = ga       (x,varargin)
    [x]   = arsa     (x,varargin)
    [x,y] = lopi     (x,y,z)
    [x,y] = sa       (x,y,varargin)
    [x,y] = fw       (x,y,z)
    %experimental/under testing
    [x,y] = test1     (x,varargin)
    [x,y] = test2     (x,y)
        
    %utilities:
    [x] = measures      (x,y)
    [x] = measures_naive(x,y) %to be removed
    [x] = template      (x,varargin)

    %permutations processing & combining
    [x]   = pdm      (x, y)
    [x]   = invperm  (x)
    [x]   = flip     (x, y)
    [x]   = matperm  (x)
    [x]   = distperm (x, y, varargin)
    [x]   = stats    (x, y, varargin)
    [x]   = relation (x)
    [x,y] = consensus(x, y, varargin)
    [x]   = cycles   (x)
    [x,y] = order    (x)
  end

end