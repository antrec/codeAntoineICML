addpath(genpath(pwd));
%%
% Generate (ground truth) similarity matrix
n = 200;
noise=0.0;
alpha=1;
% S = gen_dense_sim(n, noise, alpha);
bd = floor(n/10);
% bd = n;
nOutLim = floor(1/2 * (n + 1 - bd));
iOut = 0;
nOut = iOut*nOutLim;
iSimu=1;
rng(iSimu);
S = bandDiagOutSimMatrix(n, bd, nOut);
% S = S + 1e-0*gen_dense_sim(n,0,1);
S = S.*gen_dense_sim(n,0.,1);
% S = gen_synth_hiC_sim(n,10,0.1);
% rp = randperm(n)';
% S = S(rp,rp);
% [~,rpinv]=sort(rp);
% S = S + S';
% S = proj2RmatWithLinProg(S,1);

% (Observed) matrix with duplications
sz_ratio = 0.5;
dupl_prop = 0.5;
% bigrp = randperm(n)';

[A, Z, c] = gen_dupl_mat(S, sz_ratio, dupl_prop);
% Z = Z(bigrp,:);
A = tril(A,0) + tril(A,-1)';
% A = tril(A,-1); A = A + A';
smalln = size(A,1);
rp = randperm(smalln)';
% A = A(rp,rp);
% Z = Z(:,rp);
% c = c(rp);
figure(1); subplot(1,2,1); imagesc(S); colorbar; subplot(1,2,2); imagesc(A./(c*c')); colorbar;
% figure;imagesc(A./(c*c'));
%%
SDopts= [];

SDopts.Niter = 150;
SDopts.Ztrue = Z';
SDopts.doPlot = true;
% SDopts.doWeakR = true;
% SDopts.doQuantileNorm = true;
% SDopts.Svals = S(:);
dh = bd;
% dh = n;
SDopts.dh = dh;
thismethodopts = [];
thismethodopts.Niter = 10;
thismethodopts.Nit = 80;
thismethodopts.x_0 = (1:n)';
thismethodopts.dHuber = dh;
thismethodopts.Toeplitz = 'Huber';
thismethodopts.dH = dh;
figure;
t = clock;
rsfh = @(M) spectralEtaTrick(M, thismethodopts);
% rsfh = @(M) unconsPermOpt(M, thismethodopts);
% rsfh = @(M) seriationFAQ(M, thismethodopts);
ubval = zeros(1,n);
for idiag=1:n
    ubval(idiag) = mean(diag(S,idiag-1));
end
% ubval = 1*mean(diag(S,2));
[Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts, ubval);
et = etime(clock, t);

figure; subplot(1,2,1); imagesc(S); colorbar;subplot(1,2,2); imagesc(St);colorbar;
figure; spy(Zt(:,end:-1:1),'.'); hold on; spy(Z','rd');



%%
n=100;
opts = [];
A = gen_dense_sim(n,0.1,0.9);
ubs = zeros(1,n);
for idiag=1:n
    ubs(idiag) = mean(diag(A,idiag-1));
end
opts.ubs = ubs;
Ap = proj2Rmat2(A,opts);
figure; subplot(2,1,1); imagesc(A); colorbar; subplot(2,1,2); imagesc(Ap); colorbar;