%A Frank-Wolfe algorithm that approximately minimises the QAP problem f(P):= <P*D*P',W> starting at point P,
%where D and W are assumed symmetric. The result is a projected permutation.
%From paper: Vogelstein et al., Fast Approximate Quadratic Programming for Graph Matching, PLOSone, 2015
%D, W    : the symmetric distance and the flow matrices.
%perm    : the optimising permutation
%details : history information, such as iteration count

%JY Goulermas, 2016

function [perm, details] = fw(D, W, max_iters)

  tol = 1e-6;
  
  f    = @(X) trace(X*D*X'*W);
  n    = size(D,1);
  stop = false;
  iter = 0;
  P    = ones(n)/n; %starting point
  while ~stop
    iter  = iter + 1;
    Pprev = P;
    
    %Linear approx: h(Q):= f(P) + <g(P),Q-P>, so minimise <g(P),Q>= trace(P*D*Q'*W)
    G     = 2*W*P*D; %gradient g(P):= W'PD'+WPD
    [q,~] = lapjv(G, 1e-2); %lapjv(G) minimises trace(G*Q')=trace(A(:,q)), which is of the above form <g(P),Q>
    Q     = dma.matperm(q);

    a = f(Q-P);          %2nd order coef.
    b = trace(G'*(Q-P)); %1st order coef.
    if a > 0 %upwards
      x = min(1, max(0, -b/(2*a) ));
    else
      x = f(Q)<=f(P); %select the minimum of the two
    end

    P = (1-x)*P + x*Q;
    
    diff = norm(P-Pprev,'fro');
    stop = iter > max_iters | diff < tol;
    fprintf('\n FW algorithm: iteration= %d, with step= %1.2f, diff= %2.3f', iter, x, diff)
  end

  perm          = lapjv(-P, 1e-2);
  details.iters = iter;

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






















