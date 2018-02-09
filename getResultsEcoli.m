% Make one array of results
thisExpSumm = sprintf('EcoliResults1.txt');

fileID = fopen(thisExpSumm,'w');
fprintf(fileID,' & KT & SR & Huber & R2S & Time \\\\ \n');

A = gen_sim_from_ecoli();
nAll = floor(1/2 * nnz(A));
nInDiags = ((1:n)+1).*(2*n-(1:n))/2;
dh = find((nInDiags - nAll)>0, 1, 'first');


trunc2SUMf = @(x) truncTwoSUM(x, A, dh);

thisExpName = 'testAllMethodsEcoliIncr.m';
% Load experiment results 
load(thisExpName, 'perms', 'huberscores', 'twosumscores', 'elTimes');

% add scores to the rest
testAlgos = fieldnames(perms);
for k=1:length(testAlgos)
    myAlgo = testAlgos{k};

    thisperm = getfield(perms,myAlgo);
%     % Add dist2preR ?
%     projR = proj2RmatAll(A(thisperm,thisperm));
%     thisd2R = norm(A(thisperm,thisperm) - projR, 'fro');


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
    
    fprintf(fileID, '%s & %1.2f & %1.2f & %1.2e (%1.2e) & %1.2e (%1.2e) & %1.2e\\\\ \n',...
    myAlgo, thisKDT,...
    thisSPR, ...
    thisHuber,...
    thisTruncTwoSum,...
    thisElTime...
    );
    
end

