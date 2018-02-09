function [alph_vec, true_obj, total_viol] = sn_cd_fast(A,B,comparators,alph_vec,step,reg)
    if size(comparators,2) ~= length(alph_vec)
        fprintf('Error: number of comparators and swap_amts do not match.\n');
        return
    end

    n = size(A, 1);
    m = length(alph_vec);

    %alph_vec_rev = alph_vec(m:-1:1);
    %comparators_rev = comparators(:,m:-1:1);
    
    X = A;
    %norm(X(:))
    % loop to compute P^TAP

    % drange is the range in for us to do inverses by restoring from memory
    drad = 0.2;
    drange = [0.5 - drad 0.5 + drad];
    total_viol = numel(alph_vec( abs(alph_vec-0.5)<drad));

    col_storage = zeros(n,total_viol * 2);
    row_storage = zeros(total_viol * 2,n);

    stored = 0;

    for i = m:-1:1
        alph = alph_vec(i);
        if alph < 0.999 && alph > 0.001
            idx1 = comparators(1,i);
            idx2 = comparators(2,i);
            
            first_col = X(:,idx1);
            second_col = X(:,idx2);
            X(:,idx1) = alph * first_col + (1-alph) * second_col;
            X(:,idx2) = (1-alph) * first_col + alph * second_col;
            
            first_row = X(idx1,:);
            second_row = X(idx2,:);
            X(idx1,:) = alph * first_row + (1-alph) * second_row;
            X(idx2,:) = (1-alph) * first_row + alph * second_row;
            if alph > drange(1) && alph < drange(2)
                col_storage(:,stored+1) = first_col;
                col_storage(:,stored+2) = second_col;
                row_storage(stored+1,:) = first_row;
                row_storage(stored+2,:) = second_row;
                stored = stored + 2;
            end
        elseif alph <= 0.001
            idx1 = comparators(1,i);
            idx2 = comparators(2,i);
            
            first_col = X(:,idx1);
            X(:,idx1) = X(:,idx2);
            X(:,idx2) = first_col;
            
            first_row = X(idx1,:);
            X(idx1,:) = X(idx2,:);
            X(idx2,:) = first_row;
        end 
    end
    

    %norm(X(:))
    %X = sn_permute_row(A, comparators', alph_vec);
    %X = sn_permute_row(X', comparators', alph_vec);

    Y = B;

    %gradgrad = [];
    %grad = [];

    for i = 1:m

        alph = alph_vec(i);
        idx1 = comparators(1,i);
        idx2 = comparators(2,i);

        %% first compute M^-1 X M^-1
        if alph < 0.999 && alph > 0.001
            if alph > drange(1) && alph < drange(2)
                % in the danger zone, we just restore from memory
                X(:,idx1) = col_storage(:,stored-1);
                X(:,idx2) = col_storage(:,stored);
                X(idx1,:) = row_storage(stored-1,:);
                X(idx2,:) = row_storage(stored,:);
                stored = stored - 2;
            else
                inv_alph = alph / (2 * alph - 1);
                first_col = X(:,idx1);
                second_col = X(:,idx2);
                X(:,idx1) = inv_alph * first_col + (1-inv_alph) * second_col;
                X(:,idx2) = (1-inv_alph) * first_col + inv_alph * second_col;
                
                first_row = X(idx1,:);
                second_row = X(idx2,:);
                X(idx1,:) = inv_alph * first_row + (1-inv_alph) * second_row;
                X(idx2,:) = (1-inv_alph) * first_row + inv_alph * second_row; 
            end
        elseif alph <= 0.001
            first_col = X(:,idx1);
            X(:,idx1) = X(:,idx2);
            X(:,idx2) = first_col;
            first_row = X(idx1,:);
            X(idx1,:) = X(idx2,:);
            X(idx2,:) = first_row; 
        end

        %fprintf('%d: alph = %f, inv_alpha = %f, norm = %f\n', ...
        %        i, alph, inv_alph, norm(X(:)));

        %if (norm(X(:)) / prev_norm) > 1.02
            %[first_col second_col]
            %[X(:,idx1) X(:,idx2)]
            %[first_row ; second_row]
            %[X(idx1,:); X(idx2,:)]
        %end

        %% compute second derivative
        second_diff = (X(idx1,idx1) + X(idx2,idx2) - ...
                       X(idx1,idx2) - X(idx2,idx1)) * ...
                      (Y(idx1,idx1) + Y(idx2,idx2) - ...
                       Y(idx1,idx2) - Y(idx2,idx1)) * 2;

        %% compute first derivative
        first_diff = sum((X(:,idx1) - X(:,idx2)).* ...
                         (Y(:,idx1) - Y(:,idx2)));

        first_diff = first_diff + sum((X(idx1,:) - X(idx2,:)).* ...
                                      (Y(idx1,:) - Y(idx2,:)));

        % find ax^2 + bx + c
        %quad_term_a = second_diff / 2.0;
        %quad_term_b = first_diff - second_diff * (1-alph);
        first_diff = first_diff - second_diff * (1-alph);

        %% slap on some regularization
        %first_diff = first_diff - reg * alph;
        first_diff = first_diff - reg * (alph - 0.5);
        second_diff = second_diff - reg;
        
        %% update alph by finding exact minimizer of second derivative
        if abs(first_diff) > 0.0000001 && abs(second_diff) > 0.0000001
            % compute minimizing alph
            % the alph term is there because the function is centered there
            turning_point = alph - first_diff / second_diff;
            if second_diff < 0
                if turning_point < 0.5
                    minimizer = 1;
                else
                    minimizer = 0;
                end
            else
                minimizer = min(max(turning_point,0),1);
            end
        elseif abs(first_diff) > 0.00000001
            if first_diff < 0
                minimizer = 1;
            else
                minimizer = 0;
            end
        else
            %minimizer = alph;
            minimizer = alph;
        end

        % to correct for numeri_row errors
        %if alph == 0.5
        %    fprintf('Singularity found!\n');
        %    alph = 0.5 + realmin;
        %end

        % make sure that the new point isn't too close to 0.5
        new_alph = alph + step * (minimizer - alph);
        if abs(new_alph - 0.5) < 0.001
            %new_alph = alph + 0.9 * step * (minimizer - alph);
            if alph > 0.5
                new_alph = 0.501;
            else
                new_alph = 0.409;
            end
        end
        alph = new_alph;
        alph_vec(i) = alph;
    
        %% finally, compute MYM
        if alph < 0.999 && alph > 0.001
            first_col = Y(:,idx1);
            second_col = Y(:,idx2);
            Y(:,idx1) = alph * first_col + (1-alph) * second_col;
            Y(:,idx2) = (1-alph) * first_col + alph * second_col;
            
            first_row = Y(idx1,:);
            second_row = Y(idx2,:);
            Y(idx1,:) = alph * first_row + (1-alph) * second_row;
            Y(idx2,:) = (1-alph) * first_row + alph * second_row;
        elseif alph <= 0.001
            first_col = Y(:,idx1);
            Y(:,idx1) = Y(:,idx2);
            Y(:,idx2) = first_col;
            
            first_row = Y(idx1,:);
            Y(idx1,:) = Y(idx2,:);
            Y(idx2,:) = first_row;
        end

        % store debugging information
        %grad = [grad first_diff];
        %gradgrad = [gradgrad second_diff];
    end

    % compute objective
    %norm(X(:))
    %norm(Y(:))
    true_obj = sum(dot(X,Y));

function permuted_mat = sn_inv_permute(input_mat, compare_list, alph)
    num_comparators = size(compare_list,1);
    % flip the comparator list since we want to apply the comparators closest 
    % to the output first
    compare_list = compare_list(num_comparators:-1:1,:);
    alph = alph(num_comparators:-1:1,:);

    permuted_mat = sn_permute(input_mat, compare_list, alph);