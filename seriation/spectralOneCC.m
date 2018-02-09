function [perm, fiedler] = spectralOneCC(Ssub, doRegLap)
% Reorder first connected component of matrix Ssub (using their Fiedler
% vector)

    % separate non-connected components
    [ccs, nbComp] = conn_comp(Ssub);
    nSub = size(Ssub,1);
    if nbComp >1
        fprintf('Disconnected matrix. Performing spectral ordering on the first connected component');
        Sfirst = Ssub(ccs{1},ccs{1});
        [subperm, fiedler] = spectralOneCC(Sfirst);
        perm = (1:nSub);
        perm(ccs{1}) = ccs{1}(subperm);
        return
    end

    % Avoid trivial case where no ordering is needed
    if nSub <=2
        perm=(1:nSub)';
        fiedler = perm;
        return
    end

    % Build Laplacian matrix
    d = sum(Ssub,2);
    if issparse(Ssub)
        L2 = spdiags(d, 0, nSub, nSub) - Ssub;
        L2= L2+(1e-9)*speye(nSub);
        if nargin >1 && doRegLap
            d = 1./(sqrt(d));
            Dreg = spdiags(d, 0, nSub, nSub);
            L2 = Dreg*L2*Dreg;
        end
    else
        L2 = diag(d) - Ssub;
        L2= L2+(1e-9)*eye(nSub);
        if nargin >1 && doRegLap
            d = 1./(sqrt(d));
            Dreg = diag(d);
            L2 = Dreg*L2*Dreg;
        end
    end

    % Get second smallest eigenvector (Fiedler vector)
    if issparse(Ssub) || nSub >100
        [~,lambdamax] = eigs(L2,5,'lm');
        lambdamax = lambdamax(1,1);
        opts.tol=1.e-30;
        opts.maxit=2500;
        opts.v0=(1./sqrt(nSub))*ones(nSub,1);
        [V2,eigval] = eigs(lambdamax*speye(nSub) - L2,6,'lm',opts);

        lambdadif = eigval(2,2);
        if lambdamax - lambdadif < 1.e-8
            fprintf('non connected components! - %2.3e', lambdamax-lambdadif);
        end
        fiedler = V2(:,2);
    else
        [~, D] = eig(L2);
        dd = real(diag(D));
        lambdamax = dd(end);
        [V2,~] = eig(lambdamax*eye(nSub) - full(L2));
        fiedler = V2(:,nSub-1);
    end

    % Return the closest permutation to the Fiedler vector
    fiedler = fiedler';
    [~, perm] = sort(fiedler);
    if perm(1) > perm(end)
        perm = perm(end:-1:1);
    end
    perm = perm';
end


