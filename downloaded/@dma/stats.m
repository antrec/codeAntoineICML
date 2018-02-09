%Returns some stats for specific permutation distances.
%
%'method': the type of distance measure
%n       : permutation length (can be a vector of different sizes)

%JY Goulermas, 2014

function [S] = stats(method, n, varargin)

  if mod(length(varargin),2)
    error('dma:stats', 'missing parameters')
  end
  
  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case '???for future use'
        type = varargin{i+1};
      otherwise
        error('dma:distperm', 'unknown parameter: %s', varargin{i} )
    end
  end

  switch lower(method)

    %Original case
    case 'ppm'
      z      = n^6/15 - n^4/6 + n^2/10;
      S.max  = z/2;
      if ~mod(n,2)
        S.min = (13*n^6 - 20*n^4 + 52*n^2) / 1440  +  mod(n/2,2)/2;
      else
        warning('The correction offset is ignored here for cases n= 5,11,...')
        S.min = (13*n^6 + 10*n^4 - 23*n^2 + 000) / 1440; %000=to add offset later
      end

    case 'ppc'
      S.mean = (5*n.*(n+1)) ./ (6*(2*n.^2 - 3));
      S.var  = ( (n - 2) .* (2*n.^4 + 37*n.^3 + 42*n.^2 - 45*n - 54) ) ./ ( 18*n.*(2*n.^2 - 3).^2 .* (n - 1) );
      S.max  = 1;
      if ~mod(n,2)
        S.min  = (13*n^4 - 20*n^2 + 52 + 720*mod(n/2,2)/n^2 ) / (24*(n - 1)*(n + 1)*(2*n^2 - 3));
      else
        warning('The correction offset is ignored here for cases n= 5,11,...')
        S.min  = (13*n^2 + 23 + 000) / (24*(2*n^2 - 3)); %000=to add offset later
      end
      
      
    %Absolute case (needs proof!)
    case 'ppm-abs'
      S.max = (n^4-n^2)/6;
      if ~mod(n,2)
        S.min = (5*n^4 + 4*n^2) / 48;
      else
        S.min = (5*n^4 + 10*n^2 - 15 ) / 48;
      end
      S.max = S.max/2;
      S.min = S.min/2;

    case 'ppc-abs'
      S.max = 1;
      if ~mod(n,2)
        S.min = (5*n^2 + 4) / (8*n^2-8);
      else
        S.min = (5*n^2 + 15) / (8*n^2);
      end

    otherwise
      error('dma:stats', 'unknown distance measure: %s', type )
  end

end

%--------------------------------------------------------------------------------------------------





































