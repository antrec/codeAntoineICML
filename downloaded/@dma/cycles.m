function [out] = cycles(in)
%If input is a permutation p, it decomposes p to its cycles, otherwise it
%recomposes the given cycles to a permutation. That is: cycles(cycles(p))=p for any permutation p
%No much validation is performed!
%
%input : a row or column permutation vector, or a cell array of vectors (both, all rows or all columns)
%output: a cell array of all the cycles, or a permutation vector

%JY Goulermas, 2016

  if ~iscell(in) %permutation vector given
    out  = cell(0);
    free = true(1,length(in));
    while any(free)
      out{length(out)+1} = [];
      i                  = find(free,1);
      while free(i)
        free(i)          = false;
        out{length(out)} = [out{length(out)}, i];
        i                = in(i);
      end
    end
  else %cycles given
    n   = max(cell2mat(in));
    out = NaN(n,1);
    assert( isequal( (1:n), sort(cell2mat(in)) ) ) %cycles rows only, currently
    for c = 1 : length(in)
      out(in{c}) = circshift(in{c}(:),-1);
    end
  end
  
end