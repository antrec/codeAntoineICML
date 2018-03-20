function [simMat] = gen_synth_hiC_sim(n, nChr, noiseIntens)

simMat = zeros(n,n);
% Make equally sized chromosomes
chrSize = floor(n/nChr) + 1;
chrIdxs = (1:chrSize:n);
subR = gen_dense_sim(chrSize,0,1);
subR = subR./max(subR(:));
for iChr=1:length(chrIdxs)-1
    simMat(chrIdxs(iChr):chrIdxs(iChr+1)-1,chrIdxs(iChr):chrIdxs(iChr+1)-1) = subR;
end
if chrIdxs(end) < n
    subR = gen_dense_sim(n+1-chrIdxs(end),0,1);
    subR = subR./max(subR(:));
    simMat(chrIdxs(end):n,chrIdxs(end):n) = subR;
end

noiseMat = 0.5*rand(n);
noiseMat = noiseMat + noiseMat';
simMat = simMat + max(noiseMat(:))/max(simMat(:))*noiseIntens*noiseMat;

