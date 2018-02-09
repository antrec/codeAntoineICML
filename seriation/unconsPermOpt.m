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
MAXIT = opts.Nit;
lambda = opts.lambda;
w = opts.w_0;
w = w / norm(w);
mu = opts.mu_0;
y0 = U'*opts.x_0;
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


for it = 1 : MAXIT
    

    % Get function handle with sigmoid bias in the hyperplane containing
    % the convex envelope of the permutations
    biasedfhx = @(xvar) addSigmoidBias(fh, xvar, w, mu, lambda);
    biasedfhy = @(yvar) transfo2permhyperplane(U, yvar, biasedfhx);
    
    % optimize with minFunc
    y = minFunc(biasedfhy, yprec, minFuncOpt);

    
    % Update w and mu so that some |x_i - x_j| are in the 'l1' regime of
    % Huber
    x_t   = U*y;
    w     = x_t / norm(x_t);
    % test a tatons...
%     targetOut = 3*n;
%     Dij = abs(repmat(x_t, 1, n) - repmat(x_t', n, 1));
%     numOut = sum(Dij(:) > dHuber);
%     mu = (1 - 0.3 + 2*0.3 * (targetOut > numOut)) * mu;

    % Keep closest permutation and score
    [~, p_t] = sort(x_t);
    f_t = fh(p_t);
    fvalues(it) = f_t;
    if f_t < minVal
        minVal = f_t;
        noDecrease = 0;
    else
        noDecrease = noDecrease + 1;
    end
    
    % Re test a tatons
    df = fh((1+smalleps)*p_t) - fh((1-smalleps)*p_t);
    mu = df / ds;

    
    % Plot
    if doPlot
        subplot(2,1,1);
        imagesc(A(p_t,p_t)); title(sprintf('it %d - mu : %1.2e obj : %1.2e',it,mu, f_t))
        subplot(2,1,2); plot(x_t, 'o');
        pause(0.1);
    end
    
    yprec = U'*p_t;
    
    % Add stopping criterion
    if noDecrease > maxNoDecrease
        fprintf('Stopping after %d consecutive steps without decreasing objecive...', maxNoDecrease);
        break;
    end

end

    function options = defaultOptions(n)
        options = struct;
        options.lambda = 1 / norm(randperm(n));
        options.dHuber = inf;
        options.w_0 = randperm(n)';
        options.x_0 = randperm(n)';
        options.mu_0 = 1;
        options.Nit = 100;
        options.qtileInL2 = 0.9;
        options.doPlot = false;
        options.maxNoDecrease = 10;
        options.minFuncOpts = [];
    end

    function minFuncOpts = defaultMinFuncOptions()
        minFuncOpts = struct;
        minFuncOpts.Method = 'lbfgs';
        minFuncOpts.MaxIter = 300;
        minFuncOpts.maxFunEvals = 500;
        minFuncOpts.progTol= 1e-30;
        minFuncOpts.optTol = 1e-20;
        minFuncOpts.Display = 'iter';
    end

end