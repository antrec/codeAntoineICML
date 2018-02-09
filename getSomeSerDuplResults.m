addpath(genpath(pwd));
warning('off','all');
rng(1);
%%

fns =     {'spectral',...
    'GnCR',...
    'HuberEtaTrick',...
    'AFWTiebreakHuber',...
    'UnconsHuber',...
    'ManoptHuber'};

nSim=15;
for iSim=1:nSim
    for sz_ratio=[0.5, 0.8]
        n = floor(100/sz_ratio);
        for dupl_prop = [0.8, 0.3]
            for sparsprop = [0.2, 1]
            rng(iSim);
            % Generate (ground truth) similarity matrix
%             n = 100;
            noise=0.0;
            alpha=1;
            S = gen_dense_sim(n, noise, alpha);

            % (Observed) matrix with duplications
            sz_ratio = 0.5;
            dupl_prop = 0.8;
            [A, Z, c] = gen_dupl_mat(S, sz_ratio, dupl_prop);
            dh = n;
            fn = sprintf('exp1SerDuplSimu_%d_szr%d_dp%d_sparsprop%d.mat',iSim,10*sz_ratio,10*dupl_prop, 10*sparsprop);
            load(fn,'Zs','Ss','elTimes','Z');
            thisfileName = sprintf('exp1SerDuplSimu_%d_szr%d_dp%d_sparsprop%d.txt',iSim,10*sz_ratio,10*dupl_prop, 10*sparsprop);
            printSerDupliResults(thisfileName, Zs, Ss, elTimes, Z);

            

            bd = floor(n/10);
            dh = bd;
            nOutLim = floor(1/2 * (n + 1 - bd));
            S = bandDiagOutSimMatrix(n, bd, 0);
            % S = full(gen_diag_plus_out_sim(n, floor(n/10), 0.0));
            % (Observed) matrix with duplications
            sz_ratio = 0.5;
            dupl_prop = 0.8;
            [A, Z, c] = gen_dupl_mat(S, sz_ratio, dupl_prop);
            optSD = [];

            [Zs,Ss, elTimes] = testAllMethodsDupli(A, c, dh, SDopts);
            fn = sprintf('exp2SerDuplSimu_%d_szr%d_dp%d_sparsprop%d.mat',iSim,10*sz_ratio,10*dupl_prop, 10*sparsprop);
            save(fn,'Zs','Ss','elTimes','Z');

            S = bandDiagOutSimMatrix(n, bd, nOutLim);
            % (Observed) matrix with duplications
            sz_ratio = 0.5;
            dupl_prop = 0.8;
            [A, Z, c] = gen_dupl_mat(S, sz_ratio, dupl_prop);
            optSD = [];

            [Zs,Ss, elTimes] = testAllMethodsDupli(A, c, dh, SDopts);
            fn = sprintf('exp3SerDuplSimu_%d_szr%d_dp%d_sparsprop%d.mat',iSim,10*sz_ratio,10*dupl_prop, 10*sparsprop);
            save(fn,'Zs','Ss','elTimes','Z');

            S = bandDiagOutSimMatrix(n, bd, 2*nOutLim);
            % (Observed) matrix with duplications
            sz_ratio = 0.5;
            dupl_prop = 0.8;
            [A, Z, c] = gen_dupl_mat(S, sz_ratio, dupl_prop);
            optSD = [];

            [Zs,Ss, elTimes] = testAllMethodsDupli(A, c, dh, SDopts);
            fn = sprintf('exp4SerDuplSimuSimu_%d_szr%d_dp%d_sparsprop%d.mat',iSim,10*sz_ratio,10*dupl_prop, 10*sparsprop);
            save(fn,'Zs','Ss','elTimes','Z');
            end
        end
    end
end