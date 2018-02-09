%EXPERIMENTAL!!!!!!!!!!!!!!!!!!!! do not use

%JY Goulermas, 2016

function [perm, details] = cube(D, W)

  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function [P] = g(s)
    [~,p] = sort(s);
    P     = dma.matperm(p);
  end
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function [obj] = objective(s)
    [~,p] = sort(s);
    obj   = sum(sum(D(p,p).*W));
  end
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  max_iters = 500;
  tol       = 1e-6;
  
  f    = @(X) trace(X*D*X'*W);
  n    = size(D,1);
  stop = false;
  iter = 0;
  s    = rand(n,1); %starting point
  while ~stop
    iter  = iter + 1;
    s_prev = s;

    %Linear approx: h(Q):= f(P) + <g(P),Q-P>, so minimise <g(P),Q>= trace(P*D*Q'*W)
    S     = g(s);
    G     = 2*W*S*D; %gradient g(P):= W'PD'+WPD
    %[q,~] = lapjv(G, 1e-2); %lapjv(G) minimises trace(G*Q')=trace(A(:,q)), which is of the above form <g(P),Q>
    %Q     = dma.matperm(q);
    
    options               = optimoptions(@patternsearch);
    options.Display       = 'final';
    options.MaxIterations = 222;
    q = patternsearch(@(q)sum(sum(G.*S)), rand(n,1)+0*s, [],[],[],[], zeros(n,1), ones(n,1), options);
    Q = g(q);  

    a = f(Q-P);          %2nd order coef.
    b = trace(G'*(Q-S)); %1st order coef.
    if a > 0 %upwards
      x = min(1, max(0, -b/(2*a) ));
    else
      x = f(Q)<=f(S); %select the minimum of the two
    end

    P = (1-x)*S + x*Q;
    
    diff = norm(P-Pprev,'fro');
    stop = iter > max_iters | diff < tol;
    fprintf('\n FW algorithm: iteration= %d, with step= %1.2f, diff= %2.3f', iter, x, diff)
  end

  perm          = lapjv(-P, 1e-2);
  details.iters = iter;

  
  
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






















