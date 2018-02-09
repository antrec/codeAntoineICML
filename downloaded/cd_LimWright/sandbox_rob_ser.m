%% preamble
clear all;
close all;
randSeed=1;
rng(randSeed);
addpath('/Users/antlaplante/MVA/Stage_M2/code/');
addpath(genpath('/Users/antlaplante/THESE'));

%% Rough test of QAP algo from Lim & Wright
figure(1);
% generate random matrix W
% W = rand(n);
% W = W + W';
W = gen_dense_sim(n,0,1);
figure(1); subplot(2,2,1); imagesc(W);

% generate permuted version
ptrue = randperm(n)';
[~,invp]=sort(ptrue);
D = -W(invp,invp);
% add noise to permuted matrix D
noise_lvl = 0.;
noise_mat = rand(n);
noise_mat = noise_lvl*(noise_mat + noise_mat');
D = D + noise_mat;

iter = 500;
heur='full';
reg_limit=1e8;
%
[best_obj, perm, time, iters_taken, Dnew] = cd_process(W, D, iter, heur,reg_limit);
% [~,invperm]=sort(perm);
invperm=perm;
KTreg = corr(ptrue,invperm,'type','Kendall');
KTrev = corr(ptrue,invperm(end:-1:1),'type','Kendall');
[KendallTau,revorreg] = max(KTreg,Ktrev);
if revorreg == 2
    invperm = invperm(end:-1:1);
end
subplot(2,2,3); plot(ptrue,invperm,'kx'); title(sprintf('Kendall-Tau: %1.3f',KendallTau));
% subplot(2,2,3); plot(ptrue,invperm,'.', ptrue,invperm(end:-1:1),'kx');
% subplot(2,2,3);plot(invp,perm,'r.');
subplot(2,2,2); imagesc(-D(invperm,invperm));

%% Lazy test of Robust Seriation using out-of-the-box QAP
% close(2);
figure(2);

% generate noisy serial matrix A
n=160;
noise=0.;
alpha=0.9;
% A = gen_dense_sim(n,noise,alpha);
A = eye(n);
B = gen_diag_plus_out_sim(n, floor(n/20), 0.5);
% A = A.*B;
A = A + max(A(:))*B - A.*spones(B);
A = full(A);
% A = A + max(A(:))*full(gen_diag_plus_out_sim(n, floor(1/10), 0.5));
% A = reshape((1:n^2),n,n); A = A + A';
% A = full(A);
%
% figure;
imagesc(A);
ptrue=(1:n);
ptrue=randperm(n);
[~,invp]=sort(ptrue);
A = A(invp,invp);
% subplot(2,2,1); imagesc(A);
% subplot(2,2,2); plot(invp,perm,'o');
%%
figure(2);
% use spectral ordering for initialization
p1 = spectralOneCC(A);
% p1 = (1:n)';
Aperm = A(p1,p1);
subplot(2,2,1); imagesc(Aperm);
pcomp = invp(p1);
pcomp = pcomp';
KTreg = corr((1:n)',pcomp,'type','Kendall');
KTrev = corr((1:n)',pcomp(end:-1:1),'type','Kendall');
[KendallTau,revorreg] = max([KTreg,KTrev]);
% KTreg = corr(ptrue',p1,'type','Kendall');
% KTrev = corr(ptrue',p1(end:-1:1),'type','Kendall');
% [KendallTau,revorreg] = max([KTreg,KTrev]);
    if revorreg == 2
%         subplot(2,2,2); plot(ptrue,p1(end:-1:1),'kx');
        subplot(2,2,2); plot(pcomp(end:-1:1),'kx');
%         p1 = p1(end:-1:1);
    else
%         subplot(2,2,2); plot(ptrue,p1,'kx');
        subplot(2,2,2); plot(pcomp,'kx');
    end
%     subplot(2,2,2); plot(ptrue,p1,'kx',ptrue,p1(end:-1:1),'o');
    title(sprintf('Kendall-Tau: %1.3f',KendallTau));
% subplot(2,2,2); plot(ptrue,p1,'o', ptrue,p1(end:-1:1),'xk');
subplot(2,2,3); imagesc(Aperm);
S = proj2RmatAll(Aperm);
subplot(2,2,4); imagesc(S);
[~,invperm]=sort(p1);
Sp = S(invperm,invperm);
% subplot(2,2,4); imagesc(Sp);
% pcorr = spectralOneCC(Sp);
% subplot(2,2,2); hold on; plot(ptrue,pcorr,'xr'); hold off;
% pause;

%
%% test seriationQAP2 function
figure;
T = max(1,abs(repmat((1:n)',1,n) - repmat((1:n),n,1)));
% T = min(T.^2,(n/10)^2);
D = ones(n);
for outit=1:10
    [~,perm,App] = seriationQAP2(Aperm,50,'resort');
    D = D + 1/10*T(perm,perm);
end
D = max(D(:)) - D;
figure;imagesc(D);colorbar;
figure; imagesc(Aperm.*D);

pp = spectralOneCC(Aperm.*D);
figure;plot(pp,'o');
    
%%
out_iter = 20;
iter = 1000;
heur='ct';
reg_limit=1e2;
objseq = zeros(1,out_iter);
difsums = zeros(1,out_iter);
mydiff = 1+abs(repmat((1:n)',1,n) - repmat((1:n),n,1));
% etas = 1+abs(mydiff);
etas=ones(n);

paramgm = struct;
paramgm.maxIter = 1000;
paramgm.verbose = 1;

for kout=1:out_iter
    % to use coordinate descent from Lim & Wright
    [f, perm, time, iters_taken, Aperm] = cd_process(-S, Aperm, iter, heur,reg_limit);
%     [~,Pp]=graph_matching(S,Aperm, paramgm);
%     perm = assign(Pp',1);
%     perm = Pp'*(1:n)';
%     perm = n+1-perm;
%     [~,perm]=sort(perm);
%     perm = Pp'*(1:n)';
%     [~,perm]=sort(perm);
    % to use Birkhoff polytope based relaxation for QAP
%     Aperm = A(p1,p1);
%     x0 = eye(n);
%     [f,perm,x,iter,fs,myps]=sfw(-S,Aperm,iter,x0);
%     if perm(n)<perm(1)
%         perm = perm(n:-1:1);
%     end
    subplot(2,2,1); plot(perm,'d');
    p1 = p1(perm);
    Aperm = A(p1,p1);%.*etas(p1,p1);
    S = proj2RmatAll(Aperm);
    [~,invperm]=sort(p1);
    Sp = S(invperm,invperm);
pcomp = invp(p1);
pcomp = pcomp';
KTreg = corr((1:n)',pcomp,'type','Kendall');
KTrev = corr((1:n)',pcomp(end:-1:1),'type','Kendall');
[KendallTau,revorreg] = max([KTreg,KTrev]);
% KTreg = corr(ptrue',p1,'type','Kendall');
% KTrev = corr(ptrue',p1(end:-1:1),'type','Kendall');
% [KendallTau,revorreg] = max([KTreg,KTrev]);
    if revorreg == 2
%         subplot(2,2,2); plot(ptrue,p1(end:-1:1),'kx');
        subplot(2,2,2); plot(pcomp(end:-1:1),'kx');
%         p1 = p1(end:-1:1);
    else
%         subplot(2,2,2); plot(ptrue,p1,'kx');
        subplot(2,2,2); plot(pcomp,'kx');
    end
%     subplot(2,2,2); plot(ptrue,p1,'kx',ptrue,p1(end:-1:1),'o');
    title(sprintf('Kendall-Tau: %1.3f',KendallTau));
%     subplot(2,2,2); plot(ptrue,p1,'o',ptrue,p1(end:-1:1),'xk');
%     subplot(2,2,3); imagesc(Aperm);
%     S = proj2RmatAll(Aperm);
    subplot(2,2,3); imagesc(Aperm);
    difsum = sum(sum((Aperm - S).^2));
    difsums(kout)=difsum;
    title(sprintf('d2preR:%1.2e',difsum)); 
    subplot(2,2,4); imagesc(S);
%     pcorr = spectralOneCC(Sp);
%     subplot(2,2,2); hold on; plot(invp,pcorr,'xr'); hold off;
    title(sprintf('iter %d - diffmat %1.3e',kout, full(sum(sum((Aperm - A(ptrue,ptrue)).^2)))));
    if  sum(sum((Aperm - A(ptrue,ptrue)).^2)) == 0
        break
    end
    pause(1);
%     p1 = perm;

    mydiff = 1+abs(repmat(invperm,1,n) - repmat(invperm',n,1));
%     etas = etas - (1e-2*difsum)*mydiff;
%     etas = etas + (1./(1e-2*difsum))./mydiff;
    etas = etas + (1./(1e-2*difsum))./mydiff - (1./(1e-2*difsum))*eye(n);

%     etas = etas + (1./(1e-3*difsum))./mydiff(p1,p1);


    objseq(kout) = sum(sum(abs(Aperm - S).^1));
    fprintf('obj:%1.3e',objseq(kout));
end

%% Lazy test with jovo's Birkhoff polytope relaxation for blackbox QAP
out_iter = 10;
iter = 40;
% heur='full';
reg_limit=1e9;
objseq = zeros(1,out_iter);
x0 = eye(n);
% [~,invperm]=sort(p1);
% x0 = x0(:,invperm);
difsums = zeros(1,out_iter);
etas = ones(n);
myabsdiff = 1+abs(repmat((1:n)',1,n) - repmat((1:n),n,1));
T = max(0,abs(repmat((1:n)',1,n) - repmat((1:n),n,1)));

% Check if just QAP for Robust Seriation works
x0 = eye(n);
[~,p1s]=sfw(T,Aperm,iter,x0);
[~,p2s]=sfw(T.^2,Aperm,iter,x0);
[~,psqrt]=sfw(sqrt(T),Aperm,iter,x0);
[~,pr2]=sfw(min(T.^2,(n/10)^2),Aperm,iter,x0);
figure;
subplot(2,2,1); imagesc(Aperm(p1s,p1s));
d1s = dist2R(Aperm(p1s,p1s),'fro');
title(sprintf('QAP 1SUM -%1.4e',d1s));
subplot(2,2,2); imagesc(Aperm(p2s,p2s));
d2s = dist2R(Aperm(p2s,p2s),'fro');
title(sprintf('QAP 2SUM -%1.4e',d2s));
subplot(2,2,3); imagesc(Aperm(psqrt,psqrt));
dsqrt = dist2R(Aperm(psqrt,psqrt),'fro');
title(sprintf('QAP 1/2SUM -%1.4e',dsqrt));
subplot(2,2,4); imagesc(Aperm(pr2,pr2));
[d2rs,S] = dist2R(Aperm(pr2,pr2),'fro');
title(sprintf('QAP Rbst2S -%1.4e',d2rs));
pause;

Aperm = Aperm(pr2,pr2);
p1 = p1(pr2);
%%
figure;
for kout=1:out_iter
    % to use coordinate descent from Lim & Wright
%     [best_obj, perm, time, iters_taken, Aperm] = cd_process(-S, Aperm, iter, heur,reg_limit);
    % to use Birkhoff polytope based relaxation for QAP
%     Aperm = A(p1,p1);
    x0 = 1*eye(n) + 0.*(1/n)*ones(n);
    [f,perm,x,~,fs,myps]=sfw(-S,Aperm,iter,x0);
%     if perm(n)<perm(1)
%         perm = perm(n:-1:1);
%     end
    [Pmat,Qmat] = unstack(x,n,n);
%     x0 = Pmat*x0;
    subplot(2,2,1); plot(perm,'d');
    p1 = p1(perm);
    Aperm = A(p1,p1);
%     
%     Aperm = Aperm.*(1./myabsdiff);
%     
    
    S = proj2RmatSparse(Aperm,n/2);
    [~,invperm]=sort(p1);
    Sp = S(invperm,invperm);
pcomp = invp(p1);
pcomp = pcomp';
KTreg = corr((1:n)',pcomp,'type','Kendall');
KTrev = corr((1:n)',pcomp(end:-1:1),'type','Kendall');
[KendallTau,revorreg] = max([KTreg,KTrev]);
% KTreg = corr(ptrue',p1,'type','Kendall');
% KTrev = corr(ptrue',p1(end:-1:1),'type','Kendall');
% [KendallTau,revorreg] = max([KTreg,KTrev]);
    if revorreg == 2
%         subplot(2,2,2); plot(ptrue,p1(end:-1:1),'kx');
        subplot(2,2,2); plot(pcomp(end:-1:1),'kx');
%         p1 = p1(end:-1:1);
    else
%         subplot(2,2,2); plot(ptrue,p1,'kx');
        subplot(2,2,2); plot(pcomp,'kx');
    end
%     subplot(2,2,2); plot(ptrue,p1,'kx',ptrue,p1(end:-1:1),'o');
    title(sprintf('Kendall-Tau: %1.3f',KendallTau));
%     subplot(2,2,2); plot(ptrue,p1,'o',ptrue,p1(end:-1:1),'xk');
    subplot(2,2,3); imagesc(Aperm);
    difsum = norm(Aperm-S,'fro');
%     difsum = sum(sum((Aperm - S).^2));
    difsums(kout)=difsum;
    title(sprintf('d2preR:%1.4e',difsum)); 
    
    etas = etas + (1./(1e3*difsum))./myabsdiff(p1,p1) - (1./(1e-3*difsum))*eye(n);
%     etas = etas + (1e-3)./(difsum*myabsdiff(p1,p1));

%     S = S - 1e-2*min(T.^2,(n/10)^2);

    
%     S = proj2RmatAll(Aperm);
%     subplot(2,2,3); imagesc(Aperm);
    subplot(2,2,4); imagesc(S);
%     pcorr = spectralOneCC(Sp);
%     subplot(2,2,2); hold on; plot(invp,pcorr,'xr'); hold off;
    title(sprintf('iter %d - diffmat %1.3e',kout, full(sum(sum((Aperm - A(ptrue,ptrue)).^2)))));
    if  sum(sum((Aperm - A(ptrue,ptrue)).^2)) == 0
        break
    end
    pause(1);
%     p1 = perm;


    objseq(kout) = sum(sum(abs(Aperm - S).^1));
    fprintf('obj:%1.3e',objseq(kout));
end

%%

B = Aperm - tril(Aperm,-10);
B = tril(B,0);
B = B + tril(B,-1)';
figure;imagesc(B);
B = full(spones(B));
C = B*B';
figure;imagesc(C);colorbar;
C = sparse(C);
[ii,jj,vv] = find(C);
meanv = mean(vv);
Cw = sparse(ii,jj,vv./meanv,n,n);
figure;imagesc(Cw);colorbar;
minv=min(vv)/meanv;
Aw = minv*spones(Aperm) - minv*spones(Aperm).*spones(Cw) + spones(Aperm).*Cw;
figure;imagesc(Aw);
Cs = sum(C,1); figure;imagesc(Cs);
  %% check jovo's algorithm works ...
pp=randperm(n)';
B = A(pp,pp) + 0.1*max(A(:))*rand(n);
iter=1;
[f,perm,x,~,fs,myps]=sfw(-B,A,iter,-1);
paramgm.maxIter = 30000;
% [~,Pp]=graph_matching(B,A, paramgm);
% perm = Pp'*(1:n)'; invp=perm;
[~,invp]=sort(perm);
figure;subplot(2,1,1); imagesc(A); subplot(2,1,2); imagesc(B(invp,invp));

C = B(invp,invp);

[~,Pp]=graph_matching(-T,A, paramgm);
perm = Pp'*(1:n)'; invp=perm;
figure;subplot(2,1,1); imagesc(A); subplot(2,1,2); imagesc(C(invp,invp));

%%
heur='cont';
[best_obj, perm, time, iters_taken, Ap] = cd_process(-Aperm, A, 1000, heur,reg_limit);
%%

ntest = 15;
mytimes = zeros(1,ntest);
sizes = floor(linspace(50,900,ntest));
iter = 0;
figure;
for n=sizes
    iter = iter + 1;
    t = clock;
    A = gen_dense_sim(n,0,1);
    bd = floor(0.15*n);
    mydiff = abs(repmat((1:n)',1,n) - repmat((1:n),n,1));
    Lrs = min(mydiff.^2,bd^2);
    Lrs = Lrs - bd^2;
%     Lrs = sparse(Lrs);
    A = A - triu(A,bd) - tril(A,-bd);
%     A = sparse(A);
    prm = randperm(n);
    B = A(prm,prm);
%     [f,myp,x,iter,fs,myps] = sfw(Lrs,B,10,-1);
%     figure; plot(fs); pause;
% %     [f,myp]=sfw(Lrs,B,10,eye(n));
    [f, myp] = cd_process(Lrs, B, 1000, 'resort',100);
    et = etime(clock,t);
    fprintf('QAP of size %d in %1.2es',n,et);
    subplot(2,1,1); imagesc(A);
    subplot(2,1,2); imagesc(B(myp,myp));
    title(sprintf('QAP of size %d in %1.2es',n,et));
    pause(0.1);
    mytimes(iter) = et;

end
figure; plot(sizes,mytimes);
%% test eta-trick for seriation ?
gam = 1e-2;
n = 200;
noise=0.;
alpha=0.9;
A = gen_dense_sim(n,noise,alpha);
A = eye(n);
A = A + 3.5*max(A(:))*full(gen_diag_plus_out_sim(n, floor(n/20), 0.35));
%%
At=A;
figure(3);
subplot(1,2,1); imagesc(At);
etas = ones(n);
etasprec=etas;
myeps = 0.01;
myM = 0.8;
maxiter = 30;
twosums = zeros(1,maxiter);
onesums = zeros(1,maxiter);
onesumprec=0;
twosumprec=0;

% d = sum(A./etas,2);
% L = diag(d) - A./etas;
% L= L+(1e-9)*speye(n);
% L = full(L);

myx = ones(n,1);
myetas = ones(n);
minscore = inf;
maxout=1;
for outit=1:maxout
    etas = ones(n);
    randmat = abs(rand(n));
    randmat = randmat + randmat';
    etas = etas+randmat;
    for iter=1:maxiter

        subplot(1,3,2); imagesc(etas); colorbar;
        title(sprintf('it%d 1SUM:%2.2e 2SUM:%2.2e',iter,onesumprec,twosumprec)); 
        subplot(1,3,1); imagesc(At); colorbar;
        pause(0.01);

    %     d = sum(A./etas,2);
    %     L = diag(d) - A./etas;
    %     myx = myx - gam*L*myx;
    %     [~,myp]=sort(myx);
        [myp, myx] = spectralOneCC(A./etas);
        At = A(myp,myp);
        subplot(1,3,3); plot(myx,'o');
        if iter == 1, pause(3), end;


        mydiff = repmat(myp,1,n) - repmat(myp',n,1);
%         mydiff = repmat(myx,1,n) - repmat(myx',n,1);
        myabsdiff = abs(mydiff);
        mymed = mean(myabsdiff(:));
        mymed = floor(n/10);
%         mymed = 1/2*(quantile(myabsdiff(:),1) + quantile(myabsdiff(:),0.5));
%         mymed = 1 - 1e-1;
%         etas = etas + max(myabsdiff,mymed);
        etas = 1/maxiter * (0*etas + iter*max(myabsdiff,mymed));
%         etas = 1/2 * (etas + etasprec);
%         etasprec = max(myabsdiff,mymed);
%         etas = max(myabsdiff,mymed);
    %     etas = abs(mydiff) + 1e-7;%1e-0*mean(mydiff(:));%eye(n);
%         etas = etas + abs(mydiff).*1/n;
%         if mod(iter,30)==0
%             etas = abs(mydiff) + .5*eye(n);
%         end
    %     etas = etas/max(etas(:));
    %     etas = max(myeps,abs(mydiff));
    %     etas = abs(mydiff);
    %     etas = min(myM,max(myeps,abs(mydiff)));%eye(n);
    %     gam = 1./iter;
    %     etas = etas - gam*A.*(1 - (mydiff./etas).^2);
    %     etas = etas - diag(diag(etas)) + eye(n);
    %     etas = max(0,etas);

        twosumprec = two_SUM(myp,A);
        onesumprec = huber_log_barr(myp,A,0.1,0);
        twosums(iter) = twosumprec;
        onesums(iter) = onesumprec;
        title(sprintf('it%d 1SUM:%2.2e 2SUM:%2.2e',iter,onesumprec,twosumprec)); 
        pause(0.5);
    end
    onesums = onesums(1:iter);
    twosums = twosums(1:iter);
    figure(4); plot(onesums./mean(onesums),'r'); hold on; plot(twosums./mean(twosums),'b');hold off

    if onesums(iter) < minscore
        minscore = onesums(iter);
        myetas = etas;
    end
%     myetas = myetas + etas./maxout;
    
end
[myp, myx] = spectralOneCC(A./myetas);
At = A(myp,myp);
figure;
subplot(1,3,1); imagesc(At);
subplot(1,3,2); imagesc(myetas);
subplot(1,3,3); plot(myx,'o');

%%
At=A;
subplot(1,2,1); imagesc(At);
etas = 0.5*ones(n);
myeps = 1e-3;
myM = 0.8;
maxiter = 100;
twosums = zeros(1,maxiter);
onesums = zeros(1,maxiter);
onesumprec=0;
twosumprec=0;

myx = ones(n,1);
alph = 0.02;

for iter=1:maxiter

    subplot(1,3,2); imagesc(etas); colorbar;
    title(sprintf('it%d 1SUM:%2.2e 2SUM:%2.2e',iter,onesumprec,twosumprec)); 
    subplot(1,3,1); imagesc(At); colorbar;
    pause(0.01);

%     d = sum(A./etas,2);
%     L = diag(d) - A./etas;
%     myx = myx - gam*L*myx;
%     [~,myp]=sort(myx);
    [myp, myx] = spectralOneCC(A./etas);
    At = A(myp,myp);
    subplot(1,3,3); plot(myx,'o');

    mydiff = repmat(myp,1,n) - repmat(myp',n,1);
%     mydiff = repmat(myx,1,n) - repmat(myx',n,1);
%     etas = abs(mydiff) -etas;
    etas = max(1,abs(mydiff));
%     etas = etas+ min(etas(:))+ 1e-0;%1e-0*mean(mydiff(:));%eye(n);
%     etas = alph + (1-alph)./(max(myeps,abs(mydiff)));
%     etas = 1./etas;
%     etas = etas + abs(mydiff).*1/n;
%         if mod(iter,10)==0
%             etas = abs(mydiff) + .5*eye(n);
%         end
%     etas = etas/max(etas(:));
%     etas = max(myeps,abs(mydiff));
%     etas = abs(mydiff);
%     etas = min(myM,max(myeps,abs(mydiff)));%eye(n);
%     gam = 1./iter;
%     etas = etas - gam*A.*(1 - (mydiff./etas).^2);
%     etas = etas - diag(diag(etas)) + eye(n);
%     etas = max(0,etas);

    twosumprec = two_SUM(myp,A);
    onesumprec = huber_log_barr(myp,A,0.1,0);
    twosums(iter) = twosumprec;
    onesums(iter) = onesumprec;
    title(sprintf('it%d 1SUM:%2.2e 2SUM:%2.2e',iter,onesumprec,twosumprec)); 
    pause;
end
onesums = onesums(1:iter);
twosums = twosums(1:iter);
figure(4); plot(onesums./mean(onesums),'r'); hold on; plot(twosums./mean(twosums),'b');hold off


 %% test eta-trick on real data ?
figure(5);
A = gen_sim_from_ecoli();
n = size(A,1);
[ii,jj,vv]=find(A);
At=A;
scatshowsp(At);
%
etas = ones(length(vv),1);
maxiter = 20;
twosums = zeros(1,maxiter);
onesums = zeros(1,maxiter);
onesumprec=0;
twosumprec=0;
minscore = +inf;
maxout=1;
myetas = zeros(length(vv),maxout);
for outiter=1:maxout
    for iter=1:maxiter

    %     subplot(1,2,2); scatshowsp(etas); colorbar;
        title(sprintf('it%d 1SUM:%2e2 2SUM:%2e2',iter,onesumprec,twosumprec)); 
    %     subplot(1,2,1);
        An = sparse(ii,jj,vv./etas,n,n);
        myp = spectralOneCC(An);
        At = sparse(myp(ii),myp(jj),vv,n,n);
        subplot(1,2,1); scatshowsp(At); colorbar;
        subplot(1,2,2); plot(myp,'o');
    %     At = A(myp,myp);
        etas = max(1,abs(ii-jj));
        absdiff = 1/80*max(80,abs(ii-jj));
        etas = absdiff + iter/(iter+1) * etas;
        etas = 1/maxiter * (iter*absdiff + etas);
    %     etas = abs(repmat(myp,1,n) - repmat(myp',n,1)) + eye(n);
        twosumprec = two_SUM(myp,A);
        onesumprec = huber_log_barr(myp,A,80,0);
        if onesumprec < minscore
            minscore = onesumprec;
            bestperm = myp;
        end
        twosums(iter) = twosumprec;
        onesums(iter) = onesumprec;
        title(sprintf('it%d 1SUM:%2e2 2SUM:%2e2',iter,onesumprec,twosumprec)); 
        pause(0.1);
    end
    onesums = onesums(1:iter-1);
    twosums = twosums(1:iter-1);
    figure(6); plot(onesums./mean(onesums),'dr'); hold on; plot(twosums./mean(twosums),'xb');hold off
    
    myetas(:,outiter)=etas;
end



%% mix eta-trick and "robust seriation via QAP"

% close(2);
figure(2);

% generate noisy serial matrix A
n=150;
noise=0.;
alpha=0.9;
A = gen_dense_sim(n,noise,alpha);
A = eye(n);
A = A + max(A(:))*full(gen_diag_plus_out_sim(n, floor(n/10), 0.3));
A = gen_dense_sim(n,noise,alpha);
A = eye(n);
B = gen_diag_plus_out_sim(n, floor(n/15), 0.4);
A = A + max(A(:))*B - A.*spones(B);
A = full(A);
figure;imagesc(A);
%%
ptrue=(1:n);
ptrue=randperm(n);
[~,invp]=sort(ptrue);
A = A(invp,invp);
etas = ones(n);

figure(2);
% use spectral ordering for initialization
p1 = spectralOneCC(A);
Aperm = A(p1,p1);
subplot(2,2,1); imagesc(Aperm);
pcomp = invp(p1);
pcomp = pcomp';
KTreg = corr((1:n)',pcomp,'type','Kendall');
KTrev = corr((1:n)',pcomp(end:-1:1),'type','Kendall');
[KendallTau,revorreg] = max([KTreg,KTrev]);
% KTreg = corr(ptrue',p1,'type','Kendall');
% KTrev = corr(ptrue',p1(end:-1:1),'type','Kendall');
% [KendallTau,revorreg] = max([KTreg,KTrev]);
    if revorreg == 2
%         subplot(2,2,2); plot(ptrue,p1(end:-1:1),'kx');
        subplot(2,2,2); plot(pcomp(end:-1:1),'kx');
%         p1 = p1(end:-1:1);
    else
%         subplot(2,2,2); plot(ptrue,p1,'kx');
        subplot(2,2,2); plot(pcomp,'kx');
    end
%     subplot(2,2,2); plot(ptrue,p1,'kx',ptrue,p1(end:-1:1),'o');
    title(sprintf('Kendall-Tau: %1.3f',KendallTau));
% subplot(2,2,2); plot(ptrue,p1,'o', ptrue,p1(end:-1:1),'xk');
subplot(2,2,3); imagesc(Aperm);
S = proj2RmatAll(Aperm);
subplot(2,2,4); imagesc(S);
[~,invperm]=sort(p1);
Sp = S(invperm,invperm);
% subplot(2,2,4); imagesc(Sp);
% pcorr = spectralOneCC(Sp);
% subplot(2,2,2); hold on; plot(ptrue,pcorr,'xr'); hold off;
pause(0.1);

%%
out_iter = 30;
iter = 100;
heur='full';
reg_limit=1e9;
objseq = zeros(1,out_iter);
for kout=1:out_iter
%     Aperm = A(p1,p1);
    % divide by eta
    [best_obj, perm, time, iters_taken, Aperm] = cd_process(-S./sqrt(etas(p1,p1)), Aperm, iter, heur,reg_limit);

%     [f,perm,x,~,fs,myps]=sfw(-S./sqrt(etas(p1,p1)), Aperm./sqrt(etas(p1,p1)),iter,-1);
%     if perm(n)<perm(1)
%         perm = perm(n:-1:1);
%     end
%     [Pmat,Qmat] = unstack(x,n,n);
%     Ptot = Pmat*Ptot;
%     pivec = Pmat*(1:n)';
%     x0 = Pmat*x0;
    subplot(2,2,1); plot(perm,'d');
    p1 = p1(perm);
    Aperm = A(p1,p1);
    S = proj2RmatAll(Aperm);
    [~,invperm]=sort(p1);
    Sp = S(invperm,invperm);
    
    etas = abs(repmat(invperm,1,n) - repmat(invperm',n,1)) + 1;

        
 pcomp = invp(p1);
pcomp = pcomp';
KTreg = corr((1:n)',pcomp,'type','Kendall');
KTrev = corr((1:n)',pcomp(end:-1:1),'type','Kendall');
[KendallTau,revorreg] = max([KTreg,KTrev]);
% KTreg = corr(ptrue',p1,'type','Kendall');
% KTrev = corr(ptrue',p1(end:-1:1),'type','Kendall');
% [KendallTau,revorreg] = max([KTreg,KTrev]);
    if revorreg == 2
%         subplot(2,2,2); plot(ptrue,p1(end:-1:1),'kx');
        subplot(2,2,2); plot(pcomp(end:-1:1),'kx');
%         p1 = p1(end:-1:1);
    else
%         subplot(2,2,2); plot(ptrue,p1,'kx');
        subplot(2,2,2); plot(pcomp,'kx');
    end
%     subplot(2,2,2); plot(ptrue,p1,'kx',ptrue,p1(end:-1:1),'o');
    title(sprintf('Kendall-Tau: %1.3f',KendallTau));
%     subplot(2,2,2); plot(ptrue,p1,'o',ptrue,p1(end:-1:1),'xk');
    subplot(2,2,3); imagesc(Aperm);
    S = proj2RmatAll(Aperm);
    subplot(2,2,3); imagesc(Aperm);
    subplot(2,2,4); imagesc(S);
%     pcorr = spectralOneCC(Sp);
%     subplot(2,2,2); hold on; plot(invp,pcorr,'xr'); hold off;
    title(sprintf('iter %d - diffmat %1.3e',kout, full(sum(sum((Aperm - A(ptrue,ptrue)).^2)))));
    if  sum(sum((Aperm - A(ptrue,ptrue)).^2)) == 0
        break
    end
    
%     figure(1);
    imagesc(etas);colorbar;
    pause(1);
%     p1 = perm;


    objseq(kout) = sum(sum(abs(Aperm - S).^1));
    fprintf('obj:%1.3e',objseq(kout));
end