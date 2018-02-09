%OLD SLOWER VERSION BELOW; RETAINED FOR COMPARISON PURPOSES.....
function [perm] = bea_naive(D, mode)

  n    = size(D,1);
  perm = round(n/2); %as good as a random

  while numel(perm) < n
    remaining = setdiff(1:n, perm);
    scores    = NaN( numel(remaining), numel(perm)+1 ); %store all scores in a remaining_rows X selected_rows+1 array
    for pos = 0 : numel(perm)    %for all possible positions
      for row = 1 : numel(remaining) %for all rows not selected so far
        trial_idx         = [perm(1:pos), remaining(row), perm(pos+1:end) ];
        
        scores(row,pos+1) = moe_old( D(trial_idx,:) );
      end%row
    end%pos
    [~, idx]                 = max( scores(:) );
    [chosen_row, chosen_pos] = ind2sub(size(scores), idx); %ignore ties here and use first occurence in vec(scores)
    perm                     = [perm(1:chosen_pos-1), remaining(chosen_row), perm(chosen_pos:end) ];
  end%while   
  
 if mode == 2, perm = struct('row', perm, 'col', dma.bea(D', 1) ); end

end


function [score] = moe_old(A)
  m     = size(A,2);
  shift = [zeros(1,m); A(1:end-1,:)]; %vertical shift only
  score = sum(sum(A .* shift));  
end




