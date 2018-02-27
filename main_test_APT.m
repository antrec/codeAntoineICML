currentFolder = pwd;
addpath(genpath(currentFolder));

% ns = [100, 200, 500];
% nSimu = 50;
% iOuts = [1,2,5,10,15,20];

ns = [20];
nSimu=2;
iOuts=[2];

for n = ns
    
    bds = [floor(n/10), floor(n/20)];
    
    allbds = repmat(bds, 1, nSimu*length(iOuts));
    alliOuts = repmat(iOuts, 1, nSimu*length(bds));
    alliSimu = repmat((1:nSimu), 1, length(bds)*length(iOuts));
    
    APT_run('saveResForAllMethods',{n}, allbds, alliOuts, alliSimu);
    
end
        
