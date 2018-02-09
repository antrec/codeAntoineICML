function [Zs,Ss, elTimes] = testAllMethodsDupli(A, c, dh, SDopts)

% Initialize variables to store results
elTimes = [];
Zs = [];
Ss = [];
n = sum(c);
% Run spectral to initialize other methods
% Spectral
methodName = 'spectral';
methodNames = {};
methodOpts = {};
nMethods = 1;
methodNames{nMethods} = methodName;
t = clock;
rsfh = @(M) spectralOneCC(M);
[Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
et = etime(clock, t);
fprintf('spectral computed in %1.2es',et);
elTimes = setfield(elTimes,methodName,et);
Zs = setfield(Zs,methodName,Zt);
Ss = setfield(Ss, methodName, St);
% hubscore = huberh(perm);
% twoscore = twosumh(perm);
% huberscores = setfield(huberscores,methodName,hubscore);
% twosumscores = setfield(twosumscores,methodName,twoscore);
% subplot(2,3,1); imagesc(A((1:n)',(1:n)')); title('spectral'); pause(1);
% (1:n)' = perm;
%
% Store methods names and options in cells


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
thismethodopts.x_0 = (1:n)';
thismethodopts.S_0 = (1:n)';
thismethodopts.alpha_0 = 1;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'AFWTiebreakHuber';
thismethodopts = [];
thismethodopts.x_0 = (1:n)';
thismethodopts.S_0 = (1:n)';
thismethodopts.alpha_0 = 1;
thismethodopts.dHuber = 1e-5;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;


thismethodname = 'FAQ2SUM';
thismethodopts = [];
thismethodopts.p_0 = (1:n)';
thismethodopts.Toeplitz = '2SUM';
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'FAQHuber';
thismethodopts = [];
thismethodopts.p_0 = (1:n)';
thismethodopts.Toeplitz = 'Huber';
thismethodopts.dH = dh;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'FAQR2SUM';
thismethodopts = [];
thismethodopts.p_0 = (1:n)';
thismethodopts.Toeplitz = 'R2SUM';
thismethodopts.dH = dh;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'LWCD2SUM';
thismethodopts = [];
thismethodopts.p_0 = (1:n)';
thismethodopts.Toeplitz = '2SUM';
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'LWCDHuber';
thismethodopts = [];
thismethodopts.p_0 = (1:n)';
thismethodopts.Toeplitz = 'Huber';
thismethodopts.dH = dh;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'LWCDR2SUM';
thismethodopts = [];
thismethodopts.p_0 = (1:n)';
thismethodopts.Toeplitz = 'R2SUM';
thismethodopts.dH = dh;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;


thismethodname = 'Uncons2SUM';
thismethodopts = [];
thismethodopts.w_0 = (1/norm((1:n)'))*(1:n)';
thismethodopts.x_0 = (1/norm((1:n)'))*(1:n)';
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'UnconsHuber';
thismethodopts = [];
thismethodopts.w_0 = (1/norm((1:n)'))*(1:n)';
thismethodopts.x_0 = (1/norm((1:n)'))*(1:n)';
thismethodopts.dHuber = dh;
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'Manopt2SUM';
thismethodopts = [];
thismethodopts.x_0 = (1:n)';
nMethods = nMethods + 1;
methodNames{nMethods} = thismethodname;
methodOpts{nMethods} = thismethodopts;

thismethodname = 'ManoptHuber';
thismethodopts = [];
thismethodopts.x_0 = (1:n)';
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
rsfh = @(M) dma.gncr(-M);
[Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
et = etime(clock, t);
fprintf('spectral computed in %1.2es',et);
elTimes = setfield(elTimes,methodName,et);
Zs = setfield(Zs,methodName,Zt);
Ss = setfield(Ss, methodName, St);


methodName = 'HuberEtaTrick';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
rsfh = @(M) spectralEtaTrick(M, opts);
[Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
et = etime(clock, t);
fprintf('spectral computed in %1.2es',et);
elTimes = setfield(elTimes,methodName,et);
Zs = setfield(Zs,methodName,Zt);
Ss = setfield(Ss, methodName, St);


% methodName = 'AFWTiebreak2SUM';
% iMethod = namesToIndex(methodName);
% opts = methodOpts{iMethod};
% t = clock;
% rsfh = @(M) permAwayFrankWolfe(M, opts);
% [Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
% et = etime(clock, t);
% fprintf('spectral computed in %1.2es',et);
% elTimes = setfield(elTimes,methodName,et);
% Zs = setfield(Zs,methodName,Zt);
% Ss = setfield(Ss, methodName, St);

methodName = 'AFWTiebreakHuber';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
rsfh = @(M) permAwayFrankWolfe(M, opts);
[Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
et = etime(clock, t);
fprintf('spectral computed in %1.2es',et);
elTimes = setfield(elTimes,methodName,et);
Zs = setfield(Zs,methodName,Zt);
Ss = setfield(Ss, methodName, St);

% 
% methodName = 'FAQ2SUM';
% iMethod = namesToIndex(methodName);
% opts = methodOpts{iMethod};
% t = clock;
% rsfh = @(M) seriationFAQ(M, opts);
% [Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
% et = etime(clock, t);
% fprintf('spectral computed in %1.2es',et);
% elTimes = setfield(elTimes,methodName,et);
% Zs = setfield(Zs,methodName,Zt);
% Ss = setfield(Ss, methodName, St);
% 
% methodName = 'FAQHuber';
% iMethod = namesToIndex(methodName);
% opts = methodOpts{iMethod};
% t = clock;
% rsfh = @(M) seriationFAQ(M, opts);
% [Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
% et = etime(clock, t);
% fprintf('spectral computed in %1.2es',et);
% elTimes = setfield(elTimes,methodName,et);
% Zs = setfield(Zs,methodName,Zt);
% Ss = setfield(Ss, methodName, St);
% 
% methodName = 'FAQR2SUM';
% iMethod = namesToIndex(methodName);
% opts = methodOpts{iMethod};
% t = clock;
% rsfh = @(M) seriationFAQ(M, opts);
% [Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
% et = etime(clock, t);
% fprintf('spectral computed in %1.2es',et);
% elTimes = setfield(elTimes,methodName,et);
% Zs = setfield(Zs,methodName,Zt);
% Ss = setfield(Ss, methodName, St);
% 
% methodName = 'LWCD2SUM';
% iMethod = namesToIndex(methodName);
% opts = methodOpts{iMethod};
% t = clock;
% rsfh = @(M) seriationLWCD(M, opts);
% [Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
% et = etime(clock, t);
% fprintf('spectral computed in %1.2es',et);
% elTimes = setfield(elTimes,methodName,et);
% Zs = setfield(Zs,methodName,Zt);
% Ss = setfield(Ss, methodName, St);
% 
% methodName = 'LWCDHuber';
% iMethod = namesToIndex(methodName);
% opts = methodOpts{iMethod};
% t = clock;
% rsfh = @(M) seriationLWCD(M, opts);
% [Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
% et = etime(clock, t);
% fprintf('spectral computed in %1.2es',et);
% elTimes = setfield(elTimes,methodName,et);
% Zs = setfield(Zs,methodName,Zt);
% Ss = setfield(Ss, methodName, St);
% 
% methodName = 'LWCDR2SUM';
% iMethod = namesToIndex(methodName);
% opts = methodOpts{iMethod};
% t = clock;
% rsfh = @(M) seriationLWCD(M, opts);
% [Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
% et = etime(clock, t);
% fprintf('spectral computed in %1.2es',et);
% elTimes = setfield(elTimes,methodName,et);
% Zs = setfield(Zs,methodName,Zt);
% Ss = setfield(Ss, methodName, St);
% 
% methodName = 'Uncons2SUM';
% iMethod = namesToIndex(methodName);
% opts = methodOpts{iMethod};
% t = clock;
% rsfh = @(M) unconsPermOpt(M, opts);
% [Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
% et = etime(clock, t);
% fprintf('spectral computed in %1.2es',et);
% elTimes = setfield(elTimes,methodName,et);
% Zs = setfield(Zs,methodName,Zt);
% Ss = setfield(Ss, methodName, St);

methodName = 'UnconsHuber';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
rsfh = @(M) unconsPermOpt(M, opts);
[Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
et = etime(clock, t);
fprintf('spectral computed in %1.2es',et);
elTimes = setfield(elTimes,methodName,et);
Zs = setfield(Zs,methodName,Zt);
Ss = setfield(Ss, methodName, St);

% methodName = 'Manopt2SUM';
% iMethod = namesToIndex(methodName);
% opts = methodOpts{iMethod};
% t = clock;
% rsfh = @(M) permManOpt(M, opts);
% [Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
% et = etime(clock, t);
% fprintf('spectral computed in %1.2es',et);
% elTimes = setfield(elTimes,methodName,et);
% Zs = setfield(Zs,methodName,Zt);
% Ss = setfield(Ss, methodName, St);

methodName = 'ManoptHuber';
iMethod = namesToIndex(methodName);
opts = methodOpts{iMethod};
t = clock;
rsfh = @(M) permManOpt(M, opts);
[Zt, St] = seriationDuplialtProj(A, c, rsfh, SDopts);
et = etime(clock, t);
fprintf('spectral computed in %1.2es',et);
elTimes = setfield(elTimes,methodName,et);
Zs = setfield(Zs,methodName,Zt);
Ss = setfield(Ss, methodName, St);


end