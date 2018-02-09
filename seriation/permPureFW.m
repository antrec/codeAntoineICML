function [p_t, x_t, f_t, res] = permPureFW(cost_fun, lmo_fun, n, opts)
% Frank-Wolfe in Permutathedron

% parse options and sets undefined ones to default
opts_def = defaultOptions(n);
if nargin <= 2
    opts = opts_def;
else
    opts = build_opts(opts_def, opts);
end

alf = opts.ls_alpha;
bet = opts.ls_beta;
doPlot = opts.doPlot;
A = opts.A;

it           = 1;
minf = 0;
minx = [];
% init:

x_t         = opts.x_0;

% tracking results:

fvalues = [];
gap_values = [];

fprintf('running away-steps FW, for at most %d iterations\n', opts.Tmax);

while it <= opts.Tmax
  it = it + 1; 
    
    % Compute objective and gradient
    [f_t, grad] = cost_fun(x_t);
    
    % Frank-Wolfe corner
    s_FW = lmo_fun(grad);
    d_FW = s_FW - x_t;
    d = d_FW;
    
    % duality gap:
    gap = - d_FW' * grad;

    fvalues(it-1) = f_t;
    gap_values(it-1) = gap;

    if opts.verbose
      fprintf('it = %d -  f = %g - gap=%g\n', it, f_t, gap);
    end
    
    if gap < opts.TOL
    fprintf('end of AFW: reach small duality gap (gap=%g)\n', gap);
    break
    end 
    
    
    % Line search
    step = aLineSearch(cost_fun, x_t, d, grad, 1, bet, alf);

    % FW step:    
    x_t = x_t + step * d; 

    
    % plot
    if doPlot
%         if mod(it, 1000) == 0
%             fprintf('iter %d/%d ', k, opts.N_it);
        if mod(it, 10) == 0
%             fprintf('.');
            [~,pp]=sort(x_t);
            imagesc(A(pp,pp)); title(sprintf('it %d',it));
            pause(0.1);
        end
    end

end


    function [step, fxx] = aLineSearch(funh, x, dx, g, step, bet, alf)
        fx = funh(x);
        fxx = funh(x + step*dx);
        while (fxx > fx + alf * step * g'*dx)
            step = bet*step;
            fxx = funh(x + step*dx);
        end
    end

    function [options] = defaultOptions(n)
        options.Tmax = 1e3;
        options.TOL = 1e-5;
%         c = 1./2*(n+1)*ones(1, n);
%         x0 = 1/2 * (c + randperm(n));
        x0 = randperm(n)';
        if ~(x0(1) + 1 <= x0(n))
            x0 = x0(n:-1:1);
        end
        options.x_0 = x0;
        options.ls_alpha = 1e-4;
        options.ls_beta = 0.9;
        options.doPlot = false;
        options.A = eye(n);
        options.verbose = true;
    end

res.primal = fvalues;
res.gap = gap_values;
res.x_t = x_t;

% Return permutation
[~, p_t] = sort(x_t);

end