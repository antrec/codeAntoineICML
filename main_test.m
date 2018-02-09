currentFolder = pwd;
addpath(genpath(currentFolder));

%%
ns = [100, 200];
nSimu = 50;
iOuts = [1,2,5,10,15,20];

for n = ns
    
    bds = [floor(n/10), floor(n/20)];
    
    for bd = bds
        
        nOutLim = floor(1/2 * (n + 1 - bd));
            
        for iOut=iOuts
            
            nOut = iOut*nOutLim;

            for iSimu = 1 : nSimu
                rng(iSimu);
                A = bandDiagOutSimMatrix(n, bd, nOut);
                thisExpName = sprintf('simu_n%d_bd%d_iOut%d_kSim%d.mat', ...
                    n,bd,iOut,iSimu);

                % Chose parameter dh according to number of non-zero elements of A
                nAll = floor(1/2 * nnz(A));
                nInDiags = ((1:n)+1).*(2*n-(1:n))/2;
                dh = find((nInDiags - nAll)>0, 1, 'first');

                % Run all methods
                [perms, huberscores, twosumscores, elTimes] = testAllMethods(A, dh);

                % Save results
                save(thisExpName, 'perms', 'huberscores', 'twosumscores', 'elTimes');

            end

        end
        
    end
end
        
