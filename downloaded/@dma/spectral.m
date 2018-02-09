%From papers: Ding & He, Linearized cluster assignment via spectral ordering, ICML, 2004
%             Mavroeidis & Bingham, Enhancing the stability and efficiency..., Knowl Inf Syst, 2010
%             Barnard et al., ...
%
%W         : an nxn symmetric dissimilarity (not similarity) matrix
%'method'  : the type of method to be invoked
%'ext_perm': the external reference permutation vector for one method
%'ext_mix' : the mixture coefficient within [0,1] for the 'ext_perm'
%perm      : the permutation vector

%JY Goulermas, 2012

function [perm, details] = spectral(W, varargin)

  if mod(length(varargin),2)
    error('dma:spectral', 'missing parameters')
  end
  
  n       = size(W,1);
  method  = 'ding'; %default values
  ext_mix = 1.0;

  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'method'
        method = varargin{i+1};
      case 'ext_perm'
        ext_perm = varargin{i+1};
        if ~isvector(ext_perm) || size(ext_perm,1) ~= n
          error('dma:spectral', 'ext_perm is not a column vector or not of correct length')
        end
      case 'ext_mix'
        ext_mix = varargin{i+1};
      otherwise
        error('dma:spectral', 'unknown parameter: %s', varargin{i} )
    end
  end

  switch lower(method)

    %Uses the unnormalised graph Laplacian
    case 'barnard'
      D        = diag( sum(W,2) );
      [M,L]    = eig(D-W);                  %solve EV: (D-W)*q=l*q; (0,const) is an eigenpair here
      [~,idx]  = sort(diag(L), 'ascend');
      v        = M(:,idx(end));             %last eigv (maximise here) as W contains dissimilarities
      [~,perm] = sort(v);


    %Uses the asymmetric normalised Laplacian
    case 'ding'
      W        = max(max(W))-W;           %turn it to similarity for this method only
      D        = diag( sum(W,2) );
      [M,L]    = eig(D-W, D);             %solve GEV: (D-W)*q=l*D*q; (0,const) is an eigenpair here
      [~,idx]  = sort(diag(L), 'ascend');
      v        = M(:,idx(2));             %2nd eigv to minimise, last eigv to maximise (when W contains similarities/distances)
      [~,perm] = sort(v);


    %Uses the symmetric normalised Laplacian (same as asymmetric, but with making the eigv degree weighted; D^.5*v_asym==v_sym)
    case 'ding2'
      D        = sum(W,2) .^ 0.5;
      A        = W ./ bsxfun(@times, D, D'); %or D= sum(W,2); A= diag(D)^-.5 * W * diag(D)^-.5
      [M,L]    = eig(A);                     %solve EV: (D^-.5 *W* D^-.5)*q=l*q; (1,D^.5) is an eigenpair here
      [~,idx]  = sort(diag(L), 'ascend');
      v        = M(:,idx(1));                %first eigv (minimise here) as W contains dissimilarities, or idx(end-1) for similarities
      [~,perm] = sort(v);


    %Uses the symmetric Laplacian combined with a constructed Laplacian using external reference
    case 'mavroeidis'
      D     = sum(W,2) .^ 0.5;
      Ldata = W ./ bsxfun(@times, D, D');    %D^-.5 * W * D^-.5
      
      ref = dma.invperm(ext_perm);           %invert it to get the positional info of the external reference
      if false                               %Proposed scheme in original paper for W a similarity matrix
        v0   = D / norm(D);      
        v1   = ref' - sum(ref'.*D) / sum(D); %v1 is a) orth. to v0, b) has same ordering to ext_perm, c) of equidistant components
        v1   = v1 / norm(v1);
        Lext = v0*v0' + v1*v1'/2;            %if M=eigv(Lext), then: M(:,end)-v0 == M(:,end-1)-v1 == 0
      else                                   %Modification of the above proposed scheme, to accommodate W being a
        v1   = ref / norm(ref);              %dissimilarity matrix and send r to the 1st smallest eigenvalue instead
        Lext = - v1 * v1';
      end

      Ltot     = ext_mix * Ldata + (1-ext_mix) * Lext;
      [M,L]    = eig(Ltot);
      [~,idx]  = sort(diag(L), 'ascend');
      v        = M(:,idx(1));
      [~,perm] = sort(v);
      perm     = flipud(perm); %flip to keep close to original order (needed by the modified scheme)


    otherwise
      error('dma:spectral', 'unknown method: %s', method )
  end
  
  details.eigv = v;

end

%--------------------------------------------------------------------------------------------------






































