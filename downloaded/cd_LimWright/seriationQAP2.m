function [best_obj, perm, Ap] = seriationQAP2(A, iter, heur)


    n = size(A,1);
    Ap = A;
%     heur = 'full';
    
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
    
    % test adding some 1SUM additive term to GM
    T = max(1,abs(repmat((1:n)',1,n) - repmat((1:n),n,1)));
    mylambda = 1e-0*sum(A(:))/sum(T(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % USER DEFINED PARAMS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if continuation == true
    %     reg_limit = 100; % compute largest and smallest eigenvalues for A and S ? consider they are the same. Hence the value of the Lipschitz constant is ~lambda_max^2
        [~,eigvals] = eigs(A,1,'lm');
        reg_limit = eigvals(1,1);
    else
        reg_limit = 0;
    end
    fprintf('reg_limit=%1.2e',reg_limit);
    tol = 1e-30;
    reg_increment = reg_limit / 100.0;
    if abs(reg_increment) < 0.1
        reg_increment = reg_limit * 0.9 + 0.01;
    end
    reg = 0;
    step = 1;

    fprintf('Step Length: %f, Regularization: %f\n', step, reg);
    
    resort_limit = 100;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    alphapar = 0.5;
    alpha = alphapar*ones(m,1) + (1-alphapar)*rand(m,1);
    
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
    S = proj2RmatAll(Ap);
    subplot(1,2,1); imagesc(Ap);
    subplot(1,2,2); imagesc(S);
    pause(0.1);
    obj_prec = inf;
    myeps = 1e-2;
    cpt_in = 0;
    max_in_iter = 10;
    gampar = 0.9;
    gam = 0.;
    
    for i = 1:iter
        
        % perform coordinate descent steps
        [alpha, true_obj, viols] = sn_cd_fast(-S, Ap,comparators',alpha,step,reg);
        
        unrounded_alpha = alpha;
        reg_obj = true_obj - reg * sum((alpha - 0.5 * ones(size(alpha))).^2);

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
%                     % restart alpha ?
%                     comparators = nchoosek(1:n,2);
%                     m = size(comparators,1);
%                     alpha = alphapar*ones(m,1) + (1-alphapar)*rand(m,1);
                    break
                end
                if num_resort < resort_limit
                    num_resort = num_resort + 1;
                    fprintf('Resorting the D matrix..\n');
                    sorting_birk = sn_permute(eye(n),comparators,round(alpha));
                    Ap = sorting_birk' * Ap * sorting_birk;
                    subplot(1,2,1); imagesc(Ap);
                    inter_obj1 = sum(dot(-S,Ap));
                    S = proj2RmatAll(Ap);
                    % A.R. add a 1SUM penalty
%                     S = S - mylambda*sqrt(T);
%                     S = S - mylambda*T.^1/4;
%                     S = -T;
                    S = S - mylambda*min(T.^2,(n/10)^2);
                    subplot(1,2,2); imagesc(S); colorbar;
                    inter_obj2 = sum(dot(-S,Ap));
                    fprintf('Update of S. Obj:%1.4e->%1.4e->%1.4e\n',true_obj, inter_obj1, inter_obj2);
                    pause(0.1);
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
%     S = proj2RmatAll(Ap);
    [alphax, rounded_obj, viols] = sn_cd_fast(-S,Ap,comparators',round(alpha),0,0);
    iters_taken = i;
    fprintf('Objective from final alpha rounded: %d\n',rounded_obj);

    best_obj = rounded_obj;
    fprintf('Time taken: %f\n',time);
    fprintf('=== Best Sorting Network Objective: %d \n\n', best_obj);

    sorting_birk = sn_permute(eye(n), comparators, round(alpha));
   
%     sorting_myp = assign(sorting_birk,1);
%     perm = sorting_myp;
    Ap = sorting_birk' * Ap * sorting_birk;
    % A.R. added
    % A.R. added
    mypermmat = sorting_birk'*mypermmat;
    perm = assign(mypermmat,1);
    [~,perm]=sort(perm);

