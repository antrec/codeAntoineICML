currentFolder = pwd;
addpath(genpath(currentFolder));

%%
ns = [100, 200];
nSimu = 50;
iOuts = [1,2,5,10,15,20];

fns =     {'spectral',...
    'GnCR',...
    'HuberEtaTrick',...
    'AFWTiebreak2SUM',...
    'AFWTiebreakHuber',...
    'FAQ2SUM',...
    'FAQHuber',...
    'FAQR2SUM',...
    'LWCD2SUM',...
    'LWCDHuber',...
    'LWCDR2SUM',...
    'Uncons2SUM',...
    'UnconsHuber',...
    'Manopt2SUM',...
    'ManoptHuber'};

for n = ns
    
    bds = [floor(n/10), floor(n/20)];
    trueperm = (1:n)';
    
    for bd = bds
        
        nOutLim = floor(1/2 * (n + 1 - bd));
        
        alliOutsF = sprintf('simu_n%d_bd%d.txt', ...
                    n,bd);
        fIDiOut = fopen(alliOutsF, 'w');
        fprintf(fIDiOut,' & KT & SR & Huber & R2S & D2R & Time \\\\ \n');

        for iOut=iOuts
            
            nOut = iOut*nOutLim;
            
            % Chose parameter dh according to number of non-zero elements of A
            A = bandDiagOutSimMatrix(n, bd, nOut);
            nAll = floor(1/2 * nnz(A));
            nInDiags = ((1:n)+1).*(2*n-(1:n))/2;
            dh = find((nInDiags - nAll)>0, 1, 'first');
  
            %%
            % Make one array of results
            thisExpSumm = sprintf('simu_n%d_bd%d_iOut%d_dh%d.txt', ...
                    n,bd,iOut, dh);
                
            fileID = fopen(thisExpSumm,'w');
                
            KDT = [];
            SPR = [];
            TwoSums = [];
            Hubers = [];
            TruncTwoSums = [];
            D2R = [];
            ELT = [];
            for thisfnidx=1:length(fns)
                thisfn = fns{thisfnidx};
                KDT=setfield(KDT,thisfn,[]);
                SPR=setfield(SPR,thisfn,[]);
                TwoSums=setfield(TwoSums,thisfn,[]);
                Hubers=setfield(Hubers,thisfn,[]);
                TruncTwoSums=setfield(TruncTwoSums,thisfn,[]);
                ELT=setfield(ELT,thisfn,[]);
                D2R=setfield(D2R,thisfn,[]);
            end

fprintf(fileID,' & KT & SR & Huber & R2S & D2R & Time \\\\ \n');

            for iSimu = 1 : nSimu
                rng(iSimu);
                A = bandDiagOutSimMatrix(n, bd, nOut);
                
                
                trunc2SUMf = @(x) truncTwoSUM(x, A, bd);

                thisExpName = sprintf('simu_n%d_bd%d_iOut%d_kSim%d.mat', ...
                    n,bd,iOut,iSimu);
                % Load experiment results 
                load(thisExpName, 'perms', 'huberscores', 'twosumscores', 'elTimes');
                
                % add scores to the rest
                testAlgos = fieldnames(perms);
                for k=1:length(testAlgos)
                    myAlgo = testAlgos{k};
                    
                    thisperm = getfield(perms,myAlgo);
                    % Add dist2preR ?
                    projR = proj2RmatAll(A(thisperm,thisperm));
                    thisd2R = norm(A(thisperm,thisperm) - projR, 'fro');
                    
                    
                    KTreg = corr(trueperm,thisperm,'type','Kendall');
                    KTrev = corr(trueperm,thisperm(end:-1:1),'type','Kendall');
                    [thisKDT,revorreg] = max([KTreg,KTrev]);
                    SRreg = corr(trueperm,thisperm,'type','Spearman');
                    SRrev = corr(trueperm,thisperm(end:-1:1),'type','Spearman');
                    [thisSPR,revorreg] = max([KTreg,KTrev]);  
                    thisTruncTwoSum = trunc2SUMf(thisperm);
                    thisHuber = getfield(huberscores,myAlgo);
                    thisTwoSum = getfield(twosumscores,myAlgo);
                    thisElTime = getfield(elTimes,myAlgo);
                    
                    KDT = setfield(KDT, myAlgo, [getfield(KDT,myAlgo),thisKDT]);
                    SPR = setfield(SPR, myAlgo, [getfield(SPR,myAlgo), thisSPR]);
                    TwoSums = setfield(TwoSums, myAlgo, [getfield(TwoSums,myAlgo), thisTwoSum]);
                    Hubers = setfield(Hubers, myAlgo, [getfield(Hubers,myAlgo), thisHuber]);
                    TruncTwoSums = setfield(TruncTwoSums,myAlgo,[getfield(TruncTwoSums,myAlgo), thisTruncTwoSum]);
                    ELT = setfield(ELT,myAlgo, [getfield(ELT,myAlgo), thisElTime]);
                    D2R = setfield(D2R,myAlgo,[getfield(D2R,myAlgo), thisd2R]);
                    
                   
                    if iSimu == nSimu
                        % Print results
fprintf(fileID, '%s & %1.2f (%1.2f) & %1.2f (%1.2f) & %1.2e (%1.2e) & %1.2e (%1.2e) & %1.2e (%1.2e) & %1.2e (%1.2e)\\\\ \n',...
    myAlgo, mean(getfield(KDT,myAlgo)), std(getfield(KDT,myAlgo)),...
    mean(getfield(SPR,myAlgo)), std(getfield(SPR,myAlgo)),...
    mean(getfield(Hubers,myAlgo)), std(getfield(Hubers,myAlgo)),...
    mean(getfield(TruncTwoSums,myAlgo)), std(getfield(TruncTwoSums,myAlgo)),...
    mean(getfield(D2R,myAlgo)), std(getfield(D2R,myAlgo)),...
    mean(getfield(ELT,myAlgo)), std(getfield(ELT,myAlgo))...
    );
                    end


                        
                end                    
%                 fprintf('.');
                
            end
            
            fclose(fileID);                    


%                 % Run all methods
%                 [perms, huberscores, twosumscores, elTimes] = testAllMethods(A, dh);
% 
%                 % Save results
%                 save(thisExpName, 'perms', 'huberscores', 'twosumscores', 'elTimes');
%%
fprintf(fIDiOut, '%d ', iOut);
for k=1:length(testAlgos)
    myAlgo = testAlgos{k};
    fprintf(fIDiOut,'%s & %1.2f (%1.2f) ', myAlgo, mean(getfield(KDT,myAlgo)), std(getfield(KDT,myAlgo)));
end
fprintf(fIDiOut, '\\\\ \n');
        end
        fclose(fIDiOut);
        
    end
end
        
