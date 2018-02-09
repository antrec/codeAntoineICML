function obj = run(obj)

  obj.scores = [];

  switch ( lower(obj.algo) )
    
    case 'vat'
      [obj.perm, obj.vat_ties]         = dma.vat(obj.D, obj.dist_tol);
    case 'svat'
      [obj.perm]                       = dma.svat(obj.D, obj.svat_args{:});
    case 'revat'
      [obj.perm, obj.revat_details]    = dma.revat(obj.D, obj.revat_sample);
    case 'covat'
      [obj.perm, obj.covat_details]    = dma.covat(obj.D);
    case 'covat2'
      [obj.perm]                       = dma.covat2(obj.D);
    case 'bmst'
      [obj.perm, obj.covat_details]    = dma.bmst(obj.D);
    case 'axis'
      [obj.perm, obj.axis_details]    = dma.axis(obj.D, obj.axis_args{:});
     case 'r2e'
      [obj.perm, obj.r2e_details]      = dma.r2e(obj.D, obj.mode);
    case 'spin'
      [obj.perm, obj.spin_details]     = dma.spin(obj.D, obj.spin_args{:});
    case 'corder'
      [obj.perm]                       = dma.corder(obj.D);
    case 'pca'
      [obj.perm]                       = dma.pca(obj.D, obj.mode);
    case 'mds'
      [obj.perm, obj.mds_details]      = dma.mds(obj.D, obj.mds_args{:});
    case 'svd'
      [obj.perm]                       = dma.svd(obj.D, obj.svd_args{:}, 'mode', obj.mode);
    case 'fpca'
      [obj.perm]                       = dma.fpca(obj.D, obj.fpca_args{:}, 'mode', obj.mode);
    case 'bea'
      [obj.perm]                       = dma.bea(obj.D, obj.mode);
    case 'hc'
      [obj.perm, obj.hc_details]       = dma.hc(obj.D, obj.hc_args{:});
    case 'tsp'
      [obj.perm, obj.tsp_details]      = dma.tsp(obj.D, obj.tsp_args{:});
    case 'qap'
      [obj.perm, obj.qap_details]      = dma.qap(obj.D, obj.qap_args{:});
    case 'spectral'
      [obj.perm, obj.spectral_details] = dma.spectral(obj.D, obj.spectral_args{:});
    case 'gncr'
      [obj.perm, obj.gncr_details]     = dma.gncr(obj.D, obj.gncr_args{:});
    case 'ga'
      [obj.perm, obj.ga_details]       = dma.ga(obj.D, obj.ga_args{:});
    case 'arsa'
      [obj.perm]                       = dma.arsa(obj.D, obj.arsa_args{:});
    case 'test1'
      [obj.perm, obj.test1_details]    = dma.test1(obj.D, obj.test1_args{:});
  end


  switch ( lower(obj.algo) )

    case {'covat', 'covat2', 'bmst', 'axis'} %only two-mode analysis
      obj.Q = obj.D( obj.perm.row, obj.perm.col );
      
    case {'pca', 'svd', 'bea', 'r2e', 'fpca'} %algos with ability for dual modes
      if obj.mode == 2
        obj.Q = obj.D( obj.perm.row, obj.perm.col );
      else
        obj.Q = obj.D( obj.perm, obj.perm );
      end

    otherwise  %one-mode analysis
      obj.Q = obj.D( obj.perm, obj.perm );
      
  end%switch

end




% Note1:: Matlab mapping between permutation vectors and matrices:
% [p,~,~] = find(P');          %converts row perm P*A    to A(p,:)
% [q,~,~] = find(Q);           %converts col perm A*Q    to A(:,q)
% P       = sparse(1:n, p, 1); %converts row perm A(p,:) to P*A
% Q       = sparse(q, 1:n, 1); %converts col perm A(:,q) to A*Q

