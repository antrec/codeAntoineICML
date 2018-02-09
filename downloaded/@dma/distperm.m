%Measures the distance between two ranking vectors using different metrics.
%
%p,t     : two permutation vectors that represent b different rankings for n objects.
%'method': the type of distance measure needed, followed by either 'kendall' (default), 'hamming',
%          'kendall_norm', etc.

%JY Goulermas, 2013

function [dist] = distperm(p, t, varargin)

  if isrow(p), p = p'; end
  if isrow(t), t = t'; end

  if mod(length(varargin),2)
    error('dma:distperm', 'missing parameters')
  end
  
  n    = length(p);
  type = 'kendall'; %default values

  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'type'
        type = varargin{i+1};
      otherwise
        error('dma:distperm', 'unknown parameter: %s', varargin{i} )
    end
  end

  switch lower(type)

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Counts the number of different vector elements
    case 'hamming'
      dist = sum(p~=t);
    
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Counts the total number of pairs of discordant pairs (i,j) between p & t (unnormalised)
    %i.e., the humber of pairs with p(i)<p(j) && t(i)>t(j)
    case 'kendall'
      o    = ones(n,1);
      P    = p*o';
      T    = t*o';
      dist = sum(sum(P>P' & T<T'));

      %also1 = sum(sum( abs( dma.relation(dma.invperm(p)) - dma.relation(dma.invperm(t)) ) )) / 2;
      %also2 = 0; for i = 1 : n, for j = 1 : n, also2 = also2 + ( p(i)<p(j) && t(i)>t(j)); end, end
      %P     = sign(P-P');
      %T     = sign(T-T');
      %also3 = n*(n-1)/4 - trace(P*T')/4;

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Normalised Kendall distance [-1,+1]
    case {'kendall_norm', 'tau'}
      o     = ones(n,1);
      P     = p*o';
      T     = t*o';
      disc  = sum(sum(P>P' & T<T'));
      total = n*(n-1)/2;
      dist  = 1 - 2*disc / total;
      
      %P     = sign(P-P');
      %T     = sign(T-T');
      %also1 = trace(P*T') / sqrt( trace(P*P') * trace(T*T') ); %trace(Q*Q')=n*(n-1)
      %also2 = corr(p, t, 'type', 'kendall');

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Squared Euclidean of the difference (unnormalised)
    case 'spearman'
      dist = sum( (p-t).^2 );

      %also1  = 2*n*(n+1)*(2*n+1)/6 -2 * p'*t; %const here=2*p'*p
      %also2 = sum( (p-mean(p)) .* (t-mean(t)) )*const1 + const2
      %o     = ones(n,1);
      %P     = p*o' - o*p';
      %T     = t*o' - o*t';
      %also3 = n*(n^2-1)/6 - trace(P*T')/n;

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Normalised Spearman distance [-1,+1]
    case {'spearman_norm', 'rho'}
      dist = 1 - sum( (p-t).^2 ) * 6/n/(n^2-1);

      %o     = ones(n,1);
      %P     = p*o' - o*p';
      %T     = t*o' - o*t';
      %also1 = trace(P*T') / sqrt( trace(P*P') * trace(T*T') ); %trace(Q*Q')=n^2*(n^2-1)/6
      %also2 = corr(p, t, 'type', 'spearman');
      %also3 = sum( (p-mean(p)) .* (t-mean(t)) ) / sqrt(sum((p-mean(p)).^2) * sum((t-mean(t)).^2) );

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Cayley distance (unnormalised)
    case 'cayley'
      cycles = dma.cycles( p(dma.invperm(t)) );
      dist   = n - length(cycles); %=min #cycles to make p equal to t
      
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Normalised Cayley distance [0,+1]
    case 'cayley_norm'
      dist = dma.distperm(p, t, 'type', 'cayley') / (n-1);

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Positional proximity measure
    case 'ppm'
      P    = bsxfun(@minus, p, p').^2;
      T    = bsxfun(@minus, t, t').^2;
      dist = sum(sum(P.*T)) / 2;

      %o     = ones(n,1);
      %P     = (p*o' - o*p').^2;
      %T     = (t*o' - o*t').^2;
      %also1 = trace(P*T') / 2; 

    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Normalised positional proximity measure [lb,+1]
    case {'ppm_norm', 'ppc'}
      P    = bsxfun(@minus, p, p').^2;
      T    = bsxfun(@minus, t, t').^2;
      z    = n^6/15 - n^4/6 + n^2/10; %=trace(P*P')
      dist = sum(sum(P.*T)) / z;

      %o     = ones(n,1);
      %P     = (p*o' - o*p').^2;
      %T     = (t*o' - o*t').^2;
      %also1 = trace(P*T') / sqrt( trace(P*P') * trace(T*T') );
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %Normalised positional proximity measure [0,+1]
    case {'ppc_norm'} %for even cases only!
      ppc_min = (13*n^4 - 20*n^2 + 52 + 720*mod(n/2,2)/n^2 ) / (24*(n - 1)*(n + 1)*(2*n^2 - 3));
      dist    = ( dma.distperm(p, t, 'type', 'ppc') - ppc_min ) / ( 1 - ppc_min );

    otherwise
      error('dma:distperm', 'unknown distance measure: %s', type )
  end

end

%--------------------------------------------------------------------------------------------------





































