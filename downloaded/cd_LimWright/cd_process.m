function [best_obj, perm, time, iters_taken, D] = cd_process(W, D, iter, heur,reg_limit)
    n = size(W,1);

    fprintf('========== beginning a new QAP test ==========\n');
    fprintf('Problem Size: %d\n', n);
    fprintf('Max Iterations: %d\n', iter);
    fprintf('Heuristic: %s\n', heur);

    switch heur
    case 'full'
        resort = true;
        continuation = true;
    case 'resort'
        resort = true;
        continuation = false;
    case 'cont'
        resort = false;
        continuation = true;
    otherwise
        resort = false;
        continuation = false;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % USER DEFINED PARAMS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tol = 0.001;
    reg_increment = reg_limit / 100.0;
    if abs(reg_increment) < 0.1
        reg_increment = reg_limit * 0.9 + 0.01;
    end
    reg = 0;
    step = 1;

    fprintf('Step Length: %f, Regularization: %f\n', step, reg);
    
    resort_limit = 3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    W_trans = W';
    
    comparators = bitonic_comparators(n);
    %comparators = nchoosek(1:n,2);
    m = size(comparators,1);
    %comparators = [comparators];
    comparators = [comparators; randi(n, [m 2])];
    m = m *2;

    % flip the comparator list since we want to apply the comparators closest 
    % to the output first
    comparators = comparators(m:-1:1,:);

    %alpha = 0.5 *  ones(m,1);
%     alpha = rand(m,1);
    % !!! A.R. INITIALIZATION TO IDENTITY !!!!
%     alpha =0.*alpha + 1 * ones(m,1);
     alpha = ones(m,1);
%     alpha = 0.75 * ones(m,1);
    
    alpha_prev = 0.5 * ones(m,1);
    prev_cycle_obj = inf;
    prev_rounded_val = inf;
    prev_round_resort = false;
    unrounded_alpha = alpha;
    min_rounded_obj = inf;

    num_resort = 0;
    % A.R. added
    mypermmat = eye(n);
%     perm_old = (1:n);

    tic;
    obj = 0;
    for i = 1:iter
        [alpha, true_obj, viols] = sn_cd_fast(W,D,comparators',alpha,step,reg);
        unrounded_alpha = alpha;
        reg_obj = true_obj - reg * sum((alpha - 0.5 * ones(size(alpha))).^2);
        
        % A.R. added (sparse input problem)
        if issparse(true_obj), true_obj = full(true_obj); end

        fprintf('iter %d - viol: %d, true: %.1f, reg: %.1f \n', ...
                i, viols, true_obj, reg_obj);

        % if we are likely to not be in a local minimum, we just continue
        if ((prev_cycle_obj - reg_obj)/abs(reg_obj)) > tol || ...
            abs((prev_cycle_obj - reg_obj)/(true_obj)) > 1
            prev_cycle_obj = reg_obj;
            continue
        end

        % so... we are in a LOCAL MINIMUM?
        prev_cycle_obj = inf;   % is this next one right?

        % Let us check if we are near to a permutation
        if sum((round(alpha) - alpha).^2) < (0.1^2 * n) | ~continuation
            fprintf(' At a permutation or no continuation...\n');
            sorting_birk = sn_permute(eye(n), comparators, round(alpha));
            sorting_myp = assign(sorting_birk,1);
            if resort == true
                if num_resort > 0 && sum(round(1 - alpha)) < 1
                    break
                end
                if num_resort < resort_limit
                    num_resort = num_resort + 1;
                    fprintf('Resorting the D matrix..\n');
                    sorting_birk = sn_permute(eye(n),comparators,round(alpha));
                    D = round(sorting_birk' * D * sorting_birk);
                    comparators = nchoosek(1:n,2);
                    m = size(comparators,1);
                    alpha = ones(m,1);
                    
                    % A.R. added
                    mypermmat = sorting_birk'*mypermmat;
%                     sorting_myp = assign(sorting_birk,1);
%                     perm = sorting_myp;
%                     [~,perm]=sort(perm);
%                     perm_old = perm(perm_old);
                else
                    fprintf('Terminating. Resort limit reached\n');
                    break
                end
            else
                fprintf('Terminating. No Resort allowed.\n');
                break
            end
        end

        % if we are not near a permutation yet, let us just jack up the reg
        if continuation
            reg = reg + reg_increment;
            fprintf('increasing reg amount to %f\n', reg);

            if reg >= reg_limit
                break
            end
        end
    end
    time = toc;
    [alphax, rounded_obj, viols] = sn_cd_fast(W,D,comparators',round(alpha),0,0);
    iters_taken = i;
    % A.R. added for sparse input issues
    if issparse(rounded_obj), rounded_obj = full(rounded_obj); end
    fprintf('Objective from final alpha rounded: %d\n',rounded_obj);

    best_obj = rounded_obj;
    fprintf('Time taken: %f\n',time);
    fprintf('=== Best Sorting Network Objective: %d \n\n', best_obj);

    sorting_birk = sn_permute(eye(n), comparators, round(alpha));
   
    sorting_myp = assign(sorting_birk,1);
    perm = sorting_myp;
    D = sorting_birk' * D * sorting_birk;
    % A.R. added
    % A.R. added
    mypermmat = sorting_birk'*mypermmat;
    perm = assign(mypermmat,1);
    [~,perm]=sort(perm);