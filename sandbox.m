currentFolder = pwd;
addpath(genpath(currentFolder));
clear;
close all

%%
n = 30;
A = gen_dense_sim(n,0.9,1);
A = gen_diag_plus_out_sim(n, floor(n/10), 0.08);
A = full(A);
%%
A = gen_sim_from_ecoli();
% A = A(5600:6400, 5600:6400);
A = A(1:100,1:100);


%%
% Initialize variables to store results
elTimes = [];
perms = [];
huberscores = [];
twosumscores = [];

dh = floor(n/10);

huberh = @(x) huberSUM(x, A, dh);
twosumh = @(x) two_SUM(x, A);
%%
p0 = spectralOneCC(A);
opts=[];
rp = randperm(n)';
% B = A(rp,rp);
opts.p_0 = p0;
opts.dH = dh;
pp = seriationRobGM(A,opts);
figure; subplot(1,2,1); imagesc(A(pp,pp)); subplot(1,2,2); imagesc(A(p0,p0));
%%

paramgm = struct;
paramgm.maxIter = 50000;
paramgm.verbose = 1;
rp = randperm(n)';
B = A(rp,rp);
%%

[~,Pp]=graph_matching(max(T(:))-T,B, paramgm);
%%
% perm = assign(Pp',1);
perm = Pp*(1:n)';
% perm = n+1-perm;
% [~,perm]=sort(perm);


imagesc(B(perm,perm));



%% Run spectral to initialize other methods
% Spectral
methodName = 'spectral';
t = clock;
p0 = spectralOneCC(A);
et = etime(clock, t);
fprintf('spectral computed in %1.2es',et);
elTimes.spectral = et;
perms.spectral = p0;
hubscore = huberh(p0);
twoscore = twosumh(p0);
huberscores.spectral = hubscore;
twosumscores.spectral = twoscore;
% subplot(2,3,1); imagesc(A(p0,p0)); title('spectral'); pause(1);

%%
% Store methods names and options in cells
methodNames = {};
methodOpts = {};
nMethods = 1;
methodNames{nMethods} = 'spectral';


thismethodname = 'GnCR';
thismethodopts = [];
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'HuberEtaTrick';
thismethodopts = [];
thismethodopts.dHuber = dh;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'AFWTiebreak2SUM';
thismethodopts = [];
thismethodopts.x_0 = p0;
thismethodopts.S_0 = p0;
thismethodopts.alpha_0 = 1;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'AFWTiebreakHuber';
thismethodopts = [];
thismethodopts.x_0 = p0;
thismethodopts.S_0 = p0;
thismethodopts.alpha_0 = 1;
thismethodopts.dHuber = 1e-5;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;


thismethodname = 'FAQ2SUM';
thismethodopts = [];
thismethodopts.p_0 = p0;
thismethodopts.Toeplitz = '2SUM';
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'FAQHuber';
thismethodopts = [];
thismethodopts.p_0 = p0;
thismethodopts.Toeplitz = 'Huber';
thismethodopts.dH = floor(n/10);
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'FAQR2SUM';
thismethodopts = [];
thismethodopts.p_0 = p0;
thismethodopts.Toeplitz = 'R2SUM';
thismethodopts.dH = floor(n/10);
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'LWCD2SUM';
thismethodopts = [];
thismethodopts.p_0 = p0;
thismethodopts.Toeplitz = '2SUM';
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'LWCDHuber';
thismethodopts = [];
thismethodopts.p_0 = p0;
thismethodopts.Toeplitz = 'Huber';
thismethodopts.dH = floor(n/10);
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'LWCDR2SUM';
thismethodopts = [];
thismethodopts.p_0 = p0;
thismethodopts.Toeplitz = 'R2SUM';
thismethodopts.dH = floor(n/10);
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;


thismethodname = 'Uncons2SUM';
thismethodopts = [];
thismethodopts.w_0 = (1/norm(p0))*p0;
thismethodopts.x_0 = (1/norm(p0))*p0;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'UnconsHuber';
thismethodopts = [];
thismethodopts.w_0 = (1/norm(p0))*p0;
thismethodopts.x_0 = (1/norm(p0))*p0;
thismethodopts.dHuber = floor(n/10);
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'Manopt2SUM';
thismethodopts = [];
thismethodopts.x_0 = p0;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'ManoptHuber';
thismethodopts = [];
thismethodopts.x_0 = p0;
thismethodopts.dHuber = floor(n/10);
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

% Build a mapping : method name -> index of method
namesToIndex = containers.Map();
for iMethod = 1 : nMethods
    name = methodNames{iMethod};
    namesToIndex(name) = iMethod;
end

%%
% Run methods and store results
methodName = 'GnCR';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = dma.gncr(-A);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);


methodName = 'HuberEtaTrick';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = spectralEtaTrick(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);


methodName = 'AFWTiebreak2SUM';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = permAwayFrankWolfe(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);

methodName = 'AFWTiebreakHuber';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = permAwayFrankWolfe(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);


methodName = 'FAQ2SUM';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = seriationFAQ(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);

methodName = 'FAQHuber';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = seriationFAQ(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);

methodName = 'FAQR2SUM';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = seriationFAQ(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);

methodName = 'LWCD2SUM';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = seriationLWCD(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);

methodName = 'LWCDHuber';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = seriationLWCD(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);

methodName = 'LWCDR2SUM';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = seriationLWCD(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);

methodName = 'Uncons2SUM';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = unconsPermOpt(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);

methodName = 'UnconsHuber';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = unconsPermOpt(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);

methodName = 'Manopt2SUM';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = permManOpt(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);

methodName = 'ManoptHuber';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = permManOpt(A, opts);
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);


%% Save results


%% Plot results

nCols = floor(sqrt(nMethods));
nRows = floor(nMethods/nCols) + 1;
if ~(nCols*nRows > nMethods)
    fprintf('what?');
end
figure;
for iMethod = 1 : nMethods
    methodName = methodNames{iMethod};
    thisp = getfield(perms,methodName);
    this2sumScore = getfield(twosumscores, methodName);
    thisHuberScore = getfield(huberscores, methodName);
    thisElapsedTime = getfield(elTimes, methodName);
    subplot(nRows, nCols, iMethod);
    if issparse(A)
        spy(A(thisp,thisp));
    else
        imagesc(A(thisp,thisp));
    end
    title(sprintf('%s - f_{2S}:%1.2e - f_{Hub}:%1.2e, %1.2es', ...
        methodName, this2sumScore, thisHuberScore, thisElapsedTime));
end









%%

if 0


methodNames = {'spectral', 'GnCR', 'HuberEtaTrick', ...
    'AFWTiebreak2SUM', 'AFWTiebreakHuber', 'FAQ2SUM', 'FAQHuber', ...
    'FAQR2SUM', 'LWCD2SUM', 'LWCDHuber', 'LWCDR2SUM', ...
    'Uncons2SUM', 'UnconsHuber', 'Manopt2SUM', 'ManoptHuber'};

nMethods = length(methodNames);
methodOpts = cell(1, nMethods);
for k = 1 : nMethods
    name = methodNames(k);
    namesToIndex(name) = k;
end




figure;

% Spectral
t = clock;
p0 = spectralOneCC(A);
et = etime(clock, t);
fprintf('spectral computed in %1.2es',et);
subplot(2,3,1); imagesc(A(p0,p0)); title('spectral'); pause(1);


% GnCR
t = clock;
pgncr = dma.gncr(-A);
et = etime(clock,t);
fprintf('GnCR computed in %1.2es', et);
subplot(2,3,6); imagesc(A(pgncr,pgncr)); title('GnCR'); pause(1);


% Eta-trick L1
t = clock;
opts = struct;
[pETL1, fs] = spectralEtaTrick(A, opts);
et = etime(clock, t);
fprintf('spectral eta-trick (l1) computed in %1.2es',et);
subplot(2,3,2); imagesc(A(pETL1,pETL1)); title('eta-trick l1'); pause(1);


% Eta-trick Huber
t = clock;
opts.huberDelta = floor(n/10);
[pETH, fs] = spectralEtaTrick(A, opts);
et = etime(clock, t);
fprintf('spectral eta-trick (Huber) computed in %1.2es',et);
subplot(2,3,3); imagesc(A(pETH,pETH)); title('eta-trick Huber'); pause(1);


% Away Frank-Wolfe with tie-break on 2-SUM
t = clock;
opts = struct;
opts.x_0 = p0';
opts.S_0 = p0';
opts.alpha_0 = 1;
opts.doPlot = false;
opts.verbose = false;
opts.Tmax=500;
opts.dHuber = inf;
[pAFW2, xAFW2, fAFW2, resAFW2] = permAwayFrankWolfe(A, opts);
% fh = @(x) two_SUM(x,A);
% lmoh = @(g) LmoPermuTiebreak(g);
% [pAFW2, xAFW2, fAFW2, resAFW2] = permAFW(fh, lmoh, n, opts);
et = etime(clock, t);
fprintf('A-FW 2SUM computed in %1.2es',et);
subplot(2,3,4); imagesc(A(pAFW2,pAFW2)); title('A-FW 2SUM'); pause(1);


% A-FW with tie-break on Huber
t = clock;
opts.dHuber = 1e-4;
[pAFWH, x_t, f_t, res] = permAwayFrankWolfe(A, opts);
% dh = 1e-4;
% fh = @(x) huberSUM(x,A,dh);
% [pAFWH, x_t, f_t, res] = permAFW(fh, lmoh, n, opts);
et = etime(clock, t);
fprintf('A-FW Huber computed in %1.2es',et);
subplot(2,3,5); imagesc(A(pAFWH,pAFWH)); title('A-FW Huber'); pause(1);

% QAPs with jovo (FAQ)

[pFAQ] = seriationFAQ(A, opts);
% QAPs with Lim & Wright's c.d. 
[pLWCD] = seriationLWCD(A, opts);
% Unconstrained G.D. in Permutahedron
[pUGD] = unconsPermOpt(A, opts);
% Manopt
[pManopt] = permManOpt(A, opts);

end