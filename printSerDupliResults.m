function [] = printSerDupliResults(thisfileName, Zs, Ss, elTimes, Z)


methods = {'spectral','HuberEtaTrick', 'GnCR'};

[n,N] = size(Z);
if n >= N
    Z = Z';
    bign = n;
    n = N;
    N = bign;
end


fileID = fopen(thisfileName,'w');

for iMethod=1:length(methods)
    myAlgo = methods{iMethod};
    thisZ = getfield(Zs,myAlgo);
    [n,N] = size(thisZ);
    if n >= N
    thisZ = thisZ';
    bign = n;
    n = N;
    N = bign;
    end
    thisS = getfield(Ss, myAlgo);
    thisTime = getfield(elTimes,myAlgo);
    
    thisSR = proj2RmatAll(thisS);
    thisDist2S = norm(thisSR - thisS,'fro');
%     thisDist2S = norm(thisZ*thisS*thisZ' - A,'fro');
    transpDists = zeros(1,n);
    for i=1:n
        transpDists(i) = (1/sum(Z(i,:)))* eval_twins(full(Z(i,:)),full(thisZ(i,:)));
        
%         transpDists(i) = dtw(full(Z(i,:)),full(thisZ(i,:)));
    end
%     transpDists = transpDists;
    meanTransp = mean(transpDists);
    qt = quantile(transpDists,0.75);
    stdTransp = std(transpDists);
    
    figure; spy(Z, 'k'); hold on; spy(thisZ, 'xr');
    title(thisfileName);
    
    fprintf(fileID, ...
        '%s & %1.2e & %1.2e (%1.2e) & %1.2e\\\\ \n',...
        myAlgo, thisDist2S, meanTransp, stdTransp, qt);
end
    
fclose(fileID);

