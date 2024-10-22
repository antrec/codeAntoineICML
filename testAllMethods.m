function [perms, huberscores, twosumscores, dist2Rmats, elTimes] = testAllMethods(A, dh)

% Initialize variables to store results
elTimes = [];
perms = [];
huberscores = [];
twosumscores = [];
rtwosumscores = [];
dist2Rmats = [];

% huberh = @(x) huberSUM(x, A, dh);
% twosumh = @(x) two_SUM(x, A);

n=size(A,1);
huberh = @(p) huberSUM((1:n)',A(p,p),dh);
twosumh = @(p) two_SUM((1:n)',A(p,p));
rtwosumh = @(p) truncTwoSUM((1:n)', A(p,p), dh);
getDist2R = @(p) dist2R(A(p,p));

% Run spectral to initialize other methods
% Spectral
methodName = 'spectral';
methodNames = {};
methodOpts = {};
nMethods = 1;
methodNames{nMethods} = methodName;
t = clock;
perm = spectralOneCC(A);
et = etime(clock, t);
% fprintf('spectral computed in %1.2es',et);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);
% subplot(2,3,1); imagesc(A(p0,p0)); title('spectral'); pause(1);
p0 = perm;
[~,p0inv]=sort(p0);
%
% Store methods names and options in cells


thismethodname = 'GnCR';
thismethodopts = {'method', 'twosum'};
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'HGnCR';
thismethodopts = {'method','pshuber','delta'};
thismethodopts{length(thismethodopts)+1}= dh;
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
thismethodopts.x_0 = p0inv;
thismethodopts.S_0 = p0inv;
thismethodopts.alpha_0 = 1;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'AFWTiebreakHuber';
thismethodopts = [];
thismethodopts.x_0 = p0inv;
thismethodopts.S_0 = p0inv;
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
thismethodopts.dH = dh;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'FAQR2SUM';
thismethodopts = [];
thismethodopts.p_0 = p0;
thismethodopts.Toeplitz = 'R2SUM';
thismethodopts.dH = dh;
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
thismethodopts.dH = dh;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'LWCDR2SUM';
thismethodopts = [];
thismethodopts.p_0 = p0;
thismethodopts.Toeplitz = 'R2SUM';
thismethodopts.dH = dh;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;


thismethodname = 'Uncons2SUM';
thismethodopts = [];
% thismethodopts.w_0 = (1/norm(p0inv))*p0inv;
thismethodopts.x_0 = p0inv;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'UnconsHuber';
thismethodopts = [];
% thismethodopts.w_0 = (1/norm(p0inv))*p0inv;
thismethodopts.x_0 = p0inv;
thismethodopts.dHuber = dh;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'Manopt2SUM';
thismethodopts = [];
thismethodopts.x_0 = p0inv;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'ManoptHuber';
thismethodopts = [];
thismethodopts.x_0 = p0inv;
thismethodopts.dHuber = dh;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

% Build a mapping : method name -> index of method
namesToIndex = containers.Map();
for iMethod = 1 : nMethods
    name = methodNames{iMethod};
    namesToIndex(name) = iMethod;
end

%
% Run methods and store results
methodName = 'GnCR';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = dma.gncr(-A, opts{:});
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

methodName = 'HGnCR';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
perm = dma.gncr(-A, opts{:});
et = etime(clock, t);
elTimes = setfield(elTimes,methodName,et);
perms = setfield(perms,methodName,perm);
hubscore = huberh(perm);
twoscore = twosumh(perm);
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);


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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);


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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);

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
rtwoscore = rtwosumh(perm);
thisDist2R = getDist2R(perm);
huberscores = setfield(huberscores,methodName,hubscore);
twosumscores = setfield(twosumscores,methodName,twoscore);
rtwosumscores = setfield(rtwosumscores,methodName,rtwoscore);
dist2Rmats = setfield(dist2Rmats, methodName, thisDist2R);


end