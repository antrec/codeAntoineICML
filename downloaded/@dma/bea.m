%From paper: Mc Cormick et al., Problem decomposition and data reorganization, Operations Research, 1972
%
%Implementation of the Bond Energy Algorithm with speed up.
%X   : an nxm matrix.
%mode: 1 or 2 to specify the operating mode (mode 1 acts on rows)
%perm: the permutation vector(s)

%JY Goulermas, 2012

function [perm] = bea(X, mode)
   
  n    = size(X,1);
  perm = round(n/2); %as good as a the random case: perm = [randi(n)];
  G    = X * X';     %precalculate all inner products for speedup when using MOE local contributions

  while numel(perm) < n
    remaining = setdiff(1:n, perm);
    best_row  = NaN(numel(perm)+1, 1);      %to store the best row for each given position
    best_moe  = NaN(numel(perm)+1, 1);      %to store the moe score for the above row
    for pos = 0 : numel(perm)               %scan all possible positions
      row_scores = NaN(numel(remaining),1);
      for row = 1 : numel(remaining)        %scan all rows not selected so far
        if     pos == 0,           window = [remaining(row), perm(1)                     ];
        elseif pos == numel(perm), window = [perm(end),      remaining(row)              ];
        else                       window = [perm(pos),      remaining(row), perm(pos+1) ];
        end
        row_scores(row) = moe(window, G); %speedup by using only local window to compare winnign row for a fixed position
      end%for-row
      [~,chosen_row]  = max(row_scores(:));
      total           = [ perm(1:pos), remaining(chosen_row), perm(pos+1:end) ]; %but store the total moe to be comparable with other positions later (speedup possible here as well...)
      best_row(pos+1) = chosen_row;
      best_moe(pos+1) = moe(total, G);
    end%for-pos
    [~,chosen_pos] = max(best_moe);
    chosen_row     = best_row(chosen_pos);
    perm           = [perm(1:chosen_pos-1), remaining(chosen_row), perm(chosen_pos:end) ];
  end%while   

 if mode == 2
    perm = struct('row', perm, ...
                  'col', dma.bea(X', 1) ...
                 );
  end

end

%--------------------------------------------------------------------------------------------------

%moe: calculates the measure of effectiveness (for rows only) using the precalculated G
function [score] = moe(idx, G)

  score = 0; %faster than: score = sum(G(sub2ind(size(G),idx(1:end-1),idx(2:end))));
  for k = 1 : numel(idx)-1
    score = score + G(idx(k), idx(k+1));
  end

end

