function [done] = saveResForAllMethods(n, bd, iOut, iSimu)
    done = 0;
    
    nOutLim = floor(1/2 * (n + 1 - bd));        
    nOut = iOut*nOutLim;

    rng(iSimu);
    A = bandDiagOutSimMatrix(n, bd, nOut);
    thisExpName = sprintf('simu_n%d_bd%d_iOut%d_kSim%d.mat', ...
        n,bd,iOut,iSimu);

    % Chose parameter dh according to number of non-zero elements of A
    nAll = floor(1/2 * nnz(A));
    nInDiags = ((1:n)+1).*(2*n-(1:n))/2;
    dh = find((nInDiags - nAll)>0, 1, 'first');

    % Run all methods
    [perms, huberscores, twosumscores, dist2Rmats, elTimes] = testAllMethods(A, dh);

    % Save results
    save(thisExpName, 'perms', 'huberscores', 'twosumscores', 'dist2Rmats', 'elTimes');

    pause(0.1);
    
    done = 1;
    
end
        
