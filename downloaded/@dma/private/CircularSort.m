%Sorts angles by splitting at the largest gap. Used by some seration algorithms.
%angles: a vector of angles in(-pi,+pi]
%perm  : the permutation vector from sorting angles 

%JY Goulermas & A Kostopoulos, 2012

function [perm] = CircularSort(angles)

  [val, idx] = sort(angles);
  diffs      = diff([val; val(1)+2*pi]); %treat as ring
  [~, cut]   = max(diffs);               %largest gap val(cut+1)-val(cut), scanning from -pi --> +pi
  perm       = idx([cut+1:end, 1:cut]);  %so that angles(perm(1))-angles(perm(end)) == max(diffs)

end