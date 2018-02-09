%Creates a flow matrix W for QAP minimisation with a dissimilarity matrix D, with the purpose of seriating D.
%The sign of W is adjusted to keep it a *minimisation* problem
%n   : side length of the matrix
%type: string of the name for the type of the matrix
%k   : additional parameter for: 'banded' k>0 is the band radius
%                                'exp' k>0 adjusts squashing

%JY Goulermas & A Kostopoulos, 2012

function [W] = template(n, type, k)

  if ~exist('type', 'var')
    type = 'linear';
  end
  
  [C,R]= meshgrid(1:n);
  
  switch lower(type)

    case 'squared' %squared (inertia)
      W = -(R-C).^2;

    case 'linear'%for the least squares criterion
      W = -abs(R-C); %-ve for minimisation

    case 'linear2' %from Earle2011
      W = n - abs(R-C);
      W = W - diag(diag(W)); %added to remove constant

    case 'tsp'
      W = abs(R-C)==1;

    case 'bandedX' %nonsense template; eg, maximise for large k, minimise for small :-)
      W      = abs(R-C);
      W(W>k) = 0;
      
    case 'banded' %from Earle2011; ranges from tsp (k=1) and linear2 (k=n-1);
      W      = k+1 - abs(R-C);
      W(W<0) = 0;
      W      = W - diag(diag(W)); %added to remove constant

    case 'circular'
      t      = (n - mod(n,2)) / 2;
      W      = abs(R-C);
      W(W>t) = n-mod(n,2) - W(W>t);
      W      = -W;
      
    case 'exp' %from Tsafrir2005;
      W = -(R-C).^2;
      W = exp( W ./ (n*k) );

    case 'rank1'
      u = linspace(-1, +1, n)';
      W = u*u';

    case 6%???
      W = double('?');

    otherwise
      error('dma:template', 'unknown template: %s', type )
  end
  
 W = double(W);

end