%EXPERIMENTAL!!!!!!!!!!!!!!!!!!!! do not use
%
%D   : a nxn dissimilarity matrix
%perm: the permutation vector

%JY Goulermas, 2016

function [perm, details] = test(D)
  
  details = [];
  [perm, details] = SolveLowRankQAP(D);
  
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Use spectral information from the template W to mAXmise <PDP',W> = sum(sum(W.*D(p,p)))
function [perm, iter] = SolveLowRankQAP(D)
  n = size(D,1);
  W = -dma.template(n, 'squared');
  m = rank(W); %temporary
  
  %gather 'important' eigen-pairs: l(i) & w{i}
  [M,L]   = eig(W);
  L       = diag(L);
  [~,idx] = sort(abs(L), 'descend');
  idx     = idx(1:m);
  l       = L(idx);
  w       = num2cell(M(:,idx), 1);
%ii=03; m=1;  l=l(ii); w={w{ii}}; %`TEST-TEST-TEST with spin-sts-vector ii=2
L= @(x) diag(sum(x))-x
%m=1; l=1; w={(1:n)'}; D = L(D); %full template
%m=2; l=l(2:3); w={w{2:3}}; %TEST-TEST-TEST with 2
  
  stop = false;
  iter = 0;
  %T    = ones(n)/n;
  T    =  eye(n);
  while ~stop
    iter  = iter + 1;
    Tprev = T;
    
    fprintf('\n\n\n LowRankQAP iteration %d', iter)
    
    for i = 1 : m
      %T0 = DoubleStochastic( rand(n,n) );
      P{i} = SolveRank1QAP(T, -sign(l(i))*D, w{i} ); %eigenvalue sign changes the optimisation
    end

    %a simple way of combining the m results
    weights = abs(l) / sum(abs(l));
    T       = zeros(n); %final combined permutation
    for i = 1 : m
      T = T + weights(i) * P{i}; %sum up the: cellfun(@(x,y)x*y, num2cell(weights'), P, 'UniformOutput', false);
    end

    %project to nearest permutation matrix
    perm = lapjv(-T); %need to maximise trace(T*P'), using lapjv(A) which minimizes trace(A*P')

    %keyboard
    
    diff = norm(T-Tprev,'fro');
    stop = iter > 60 | diff < 1e-13;
  end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Approximately minimise f(P):= <P*D*P',w*w'> using Frank-Wolfe, from point P.
%D is assumed symmetric. w is arbitrary.
%The matrix result can be a doubly stochastic or a permutation matrix.
function [P] = SolveRank1QAP(P, D, w)
  stop = false;
  iter = 0;
  while ~stop
    iter  = iter + 1;
    Pprev = P;

    %Linear approx: h(Q):= f(P) + <g(P),Q-P>, so minimise <g(P),Q>= w'*P*D*Q'*w
    G     = 2*w*w'*P*D; %gradient g(P):= W'PD'+WPD
    [~,r] = sort(w,      'ascend');  %Let Q=R'*S, the minimiser of w'*P*D*S'*R*w = <R*w, S*D*P'*w>
    [~,s] = sort(D*P'*w, 'descend'); %so that, R and S sort w and D*P'*w opposite directions
    q     = s(dma.invperm(r));
    Q     = dma.matperm(q); %duplicated for clarity!
    if false %some sanity-check lines
       assert( isequal(dma.matperm(r)'*dma.matperm(s), Q ) ) %Q=R'*S  
      [q_lap,~] = lapjv(G, 1e-2); %lapjv(G) minimises trace(G*Q')=trace(A(:,q)), which is of the above form <g(P),Q>
      assert( isequal(q,q_lap') )
    end

    if true
      %step search, to minimise f((1-x)*P + x*Q), for step x with 0<=x<=1.
      a = w'*(Q-P)*D*(Q-P)'*w; %= 2nd order coef. =f(Q-P)
      b = 2*w'*P*D*(Q-P)'*w;   %= 1st order coef. =trace(G'*(Q-P))
      %xx=0:.1:1; clf, plot(xx,a*xx.^2+b*xx), axis square, pause, imagesc(Q*D*Q'), pause %debug
      if a > 0 %opens upward
        x = min(1, max(0, -b/(2*a) ));
      else
        x = w'*Q*D*Q'*w <= w'*P*D*P'*w; %select the minimum from: f(Q) & f(P)
      end
    else
      %like SPIN
      x = 1;
    end

    P = (1-x)*P + x*Q;
    
    %diff = norm(P-Pprev,'fro');
    diff = norm( (P-Pprev)*D*Pprev'*w, inf);
    stop = iter > 500 | diff < 1e-7;
    fprintf('\n    Rank1QAP: iteration %d, with alpha=%1.2f, diff=%2.3f', iter, x, diff)
  end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~















































