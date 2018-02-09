function [best_obj, perm, Ap] = seriationQAPsimple(A, iter, reg_limit)


    n = size(A,1);
    Ap = A;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % USER DEFINED PARAMS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    continuation = false;
%     reg_limit = 100; % compute largest and smallest eigenvalues for A and S ? consider they are the same. Hence the value of the Lipschitz constant is ~lambda_max^2
    [~,eigvals] = eigs(A,1,'lm');
    reg_limit = eigvals(1,1)
    tol = 1e-3;
    reg_increment = reg_limit / 100.0;
    if abs(reg_increment) < 0.1
        reg_increment = reg_limit * 0.9 + 0.01;
    end
    reg = 0;
    step = 1;

    fprintf('Step Length: %f, Regularization: %f\n', step, reg);
    
    resort_limit = 10;
    
    
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
    alpha = rand(m,1);
    
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
    max_in_iter = 5;
    gampar = 0.9;
    gam = 0.;
    
    for i = 1:iter
        
        % perform coordinate descent steps
        [alpha, true_obj, viols] = sn_cd_fast(-S, Ap,comparators',alpha,step,reg);
        
        if (~((obj_prec - true_obj)/true_obj > myeps) && (cpt_in < max_in_iter) )
            obj_prec = true_obj;
            fprintf('true obj:%1.4e\n',true_obj);
            cpt_in = cpt_in +1;

            continue
        end
        
        cpt_in=0;
        
%         
%         sorting_birk = sn_permute(eye(n), comparators, round(alpha));
        sorting_birk = sn_permute(eye(n), comparators, alpha);
%         break;
%         sorting_birk = getPermMatrixFromDBQAP(sorting_birk,Ap,S,100);

        gam = gam*gampar;
        comparators = nchoosek(1:n,2);
        m = size(comparators,1);
        alpha = gam*rand(m,1) + (1-gam)*ones(m,1);

        Ap = sorting_birk' * Ap * sorting_birk;
        subplot(1,2,1); imagesc(Ap);
        
        inter_obj1 = sum(dot(-S,Ap));
        
        S = proj2RmatAll(Ap);
        subplot(1,2,2); imagesc(S);
        
        inter_obj2 = sum(dot(-S,Ap));
        
        fprintf('Update of S. Obj:%1.4e->%1.4e->%1.4e\n',true_obj, inter_obj1, inter_obj2);
        obj_prec = inter_obj2;
        
         if continuation
            reg = reg + reg_increment;
            fprintf('increasing reg amount to %f\n', reg);

            if reg >= reg_limit
                break
            end
         end
        
        pause(1);
%         fprintf('Update of S. True obj:%1.2e\n',true_obj);
        
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
    
%     getPermMatrixFromDBQAP(sorting_raw,Ap,S,nbTrial)
   
%     sorting_myp = assign(sorting_birk,1);
%     perm = sorting_myp;
    Ap = sorting_birk' * Ap * sorting_birk;
    % A.R. added
    % A.R. added
    mypermmat = sorting_birk'*mypermmat;
    perm = assign(mypermmat,1);
    [~,perm]=sort(perm);

