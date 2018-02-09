function [o,C] = order(p)
%Returns the order of a permutation
%
%p: the permutation
%o: the order, that is: dma.matperm(p)^o == eye(length(p))
%C: the cycles, in case the caller needs it

%JY Goulermas, 2016

  C = dma.cycles(p);
  o = lcm2( cellfun(@length, C) );

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%A not very efficient and maxflint-afflicted way of calculating the LCM for a set of numbers
%stored in a vector v, using the built-in lcm() for use between two integers.
function [m] = lcm2(v)
  if isscalar(v)
    m = v;
  else
    m = lcm( lcm2(v(2:end)), v(1) );
  end
end

