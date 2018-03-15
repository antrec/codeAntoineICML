function [p_t, fvalues] = unconsPermOpt(A, opts)

n = size(A,1);
opts_def = defaultOptions(n);
opts = build_opts(opts_def, opts);
minFuncOpt_def = defaultMinFuncOptions();
minFuncOpt = build_opts(minFuncOpt_def, opts.minFuncOpts);

% Build basis for hyperplane containing permutations (orthogonal to
% ones(n,1))
U = orthoConstBasis(n);
U = U(:,1:n-1);

% Init.
c = 1./2*(n+1)*ones(n,1);

MAXIT = opts.Nit;
lambda = opts.lambda;
% w = opts.w_0;
% w = w / norm(w);
% mu = opts.mu_0;
y0 = U'*(opts.x_0-c);

% Actually w0 and x0 should be the same ?
w = opts.x_0;
w = w / norm(w);
% y0 = y0 / norm(y0);
yprec = y0;
dHuber = opts.dHuber;
if dHuber == +inf
    fh = @(x) two_SUM(x, A);
    dHuber = 1;
else
    fh = @(x) huberSUM(x, A, dHuber);
end
doPlot = opts.doPlot;
fvalues = zeros(1, MAXIT);
minVal = inf;
noDecrease = 0;
maxNoDecrease = opts.maxNoDecrease;

% Constant for updating mu
smalleps = 1e-5;
ds = (exp(smalleps) - exp(-smalleps)) * exp(-1) / ((1+exp(-1 -smalleps))*(1+exp(-1+smalleps)));
df = fh((1+smalleps)/lambda*w) - fh((1-smalleps)/lambda*w);
mu = df / ds;
mu = 0.5 * mu;
% mu = 1e-3*mu;
% mu = 1e6;
% mu = 1e6;

for it = 1 : MAXIT
    

    % Get function handle with sigmoid bias in the hyperplane containing
    % the convex envelope of the permutations
    biasedfhx = @(xvar) addSigmoidBias(fh, xvar, w, mu, lambda);
    % A.R. : should remove lambda dependency
    biasedfhy = @(yvar) transfo2permhyperplane2(U, yvar, biasedfhx);
    
    % optimize with minFunc
    y = minFunc(biasedfhy, yprec, minFuncOpt);

    
    % Update w and mu so that some |x_i - x_j| are in the 'l1' regime of
    % Huber
    x_t   = U*y + c;
%     x_t = U*y*norm(U'*randperm(n)') + c;
%     w   = x_t / norm(x_t);
%     w = x_t - mean(x_t);
%     w = w / norm(x_t);
%     w = w + mean(x_t)/norm(x_t);
    % test a tatons...
%     targetOut = 3*n;
%     Dij = abs(repmat(x_t, 1, n) - repmat(x_t', n, 1));
%     Dij = Dij.*A;
% %     subplot(2,1,2); imagesc(Dij.*A);
%     numOut = sum(Dij(:) > dHuber);
%     mu = (1 - 0.3 + 2*0.3 * (targetOut > numOut)) * mu;

    % Keep closest permutation and score
    [~, p_t] = sort(x_t);
    [~, p_t] = sort(p_t);
    
%     Dij = abs(repmat(p_t, 1, n) - repmat(p_t', n, 1));
%     Dij = Dij.*A;
% %     subplot(2,1,2); imagesc(Dij.*A);
%     numOut = sum(Dij(:) > dHuber);
%     normt=norm(x_t)/norm(p_t);

    f_t = fh(p_t);
    fvalues(it) = f_t;
    if f_t < minVal
        minVal = f_t;
        minPerm = p_t;
        noDecrease = 0;
%         mu = 1.2*mu;
    else
        noDecrease = noDecrease + 1;
%         mu = mu/1.2;
    end
    
    % Re test a tatons
%     if mod(it, 1)==0
    df = fh((1+smalleps)*p_t) - fh((1-smalleps)*p_t);
%     df = fh((1+smalleps)*x_t) - fh((1-smalleps)*x_t);

    mu = 0.7*df / ds;
%     end
%     mu = 0.5 * mu;

%     mu = 1;
% if full(numOut) < 225
%     mu = mu*1.01;
% else
%     mu = mu / 1.01;
% end
% mu = 2e6;
%     mu = mu * (3*it/MAXIT);
%      mu=1e1;
% mu = 1e6;
% if it < 50
% mu = 1e-6 * (1/0.85)^it;
% % end
% mu = 1e-1;
% w = p_t / norm(p_t);
    
    % Plot
    if doPlot
        subplot(2,1,1);
        [~,pinv]=sort(p_t);
        imagesc(A(pinv,pinv)); title(sprintf('it %d - mu : %1.2e obj : %1.2e',it,mu, f_t))
        subplot(2,1,2); plot(x_t, 'o');
        pause(0.1);
    end
    
    yprec = U'*(p_t-c);
%     yprec = y;
    w = p_t / norm(p_t);
%     yprec = 1/norm(yprec) * yprec;
    
    % Add stopping criterion
    if noDecrease > maxNoDecrease
        fprintf('Stopping after %d consecutive steps without decreasing objective...', maxNoDecrease);
        break;
    end

end

[~, p_t] = sort(minPerm);


    function options = defaultOptions(n)
        options = struct;
        options.lambda = 1 / norm(randperm(n));
        options.dHuber = inf;
%         options.w_0 = randperm(n)';
        options.x_0 = randperm(n)';
        options.w_0 = options.x_0 / norm(options.x_0);
        options.mu_0 = 1;
        options.Nit = 200;
        options.qtileInL2 = 0.9;
        options.doPlot = false;
        options.maxNoDecrease = 15;
        options.minFuncOpts = [];
    end

    function minFuncOpts = defaultMinFuncOptions()
        minFuncOpts = struct;
        minFuncOpts.Method = 'lbfgs';
        minFuncOpts.MaxIter = 500;
        minFuncOpts.maxFunEvals = 1000;
        minFuncOpts.progTol= 1e-30;
        minFuncOpts.optTol = 1e-20;
        minFuncOpts.Display = 'iter';
    end

end