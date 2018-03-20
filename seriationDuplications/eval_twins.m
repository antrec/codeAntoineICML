function dist = eval_twins( z_1, z_2 )
%EVAL_TWINS Summary of this function goes here
%   Detailed explanation goes here

ind_1 = find(z_1);
ind_2 = find(z_2);

D = pdist2(ind_1, ind_2, 'cityblock');

[~, dist] = assignmentoptimal_mex(D);

end

