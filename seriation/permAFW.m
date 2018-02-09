function [p_t, x_t, f_t, res] = permAFW(cost_fun, lmo_fun, n, opts)
% Frank-Wolfe in Permutathedron

% parse options and sets undefined ones to default
opts_def = defaultOptions(n);
if nargin <= 3
    opts = opts_def;
else
    opts = build_opts(opts_def, opts);
end

alf = opts.ls_alpha;
bet = opts.ls_beta;
doPlot = opts.doPlot;
A = opts.A;

it           = 1;
minf = inf;
minx = [];
% init:

x_t         = opts.x_0;
if ~(x_t(1) + 1 <= x_t(n))
    x_t = x_t(n:-1:1);
end
S_t         = opts.S_0;
alpha_t     = opts.alpha_0;

% each column of S_t is a potential vertex
% I_active -> contains index of active vertices (same as alpha_t > 0)
mapping = containers.Map();
% this will map vertex hash to id in S_t (to see whether already seen vertex)
% alpha_t(i) == 0 implies vertex is not active anymore...
% alpha_t will be a weight vector so that x_t = S_t * alpha_t

% constructing mapping:
max_index = size(S_t,2); % keep track of size of S_t
for index = 1:max_index
    thiskey = hashing(S_t(:,index));
    mapping(thiskey) = index;
%     mapping(hashing(S_t(:,index))) = index; 
end
I_active = find(alpha_t > 0);

% tracking results:

fvalues = [];
fPermVals = [];
gap_values = [];
number_away = 0;
number_drop = 0; % counting drop steps (max stepsize for away step)


fprintf('running away-steps FW, for at most %d iterations\n', opts.Tmax);

while it <= opts.Tmax
  it = it + 1; 
    
    % Compute objective and gradient
    [f_t, grad] = cost_fun(x_t);
    
    % Frank-Wolfe corner
    s_FW = lmo_fun(grad);
    d_FW = s_FW - x_t;
    
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
    
  % away direction search:

    if ~isempty(S_t)
    id_A   = away_step(grad, S_t, I_active);
    v_A    = S_t(:, id_A);
    d_A    = x_t - v_A;
    alpha_max = alpha_t(id_A);
    else
    fprintf('error: empty support set at step (it=%f)\n', it);
    end

    % construct direction (between towards and away):

    if isempty(S_t) || - gap <= d_A' * grad
    is_aw = false; % this is a pure FW step
    d = d_FW; 
    max_step = 1;
    else
    is_aw = true; % this is an away step
    number_away = number_away+1;
    d = d_A;
    max_step = alpha_max / (1 - alpha_max);
    end
    
    % Line search
    step = aLineSearch(cost_fun, x_t, d, grad, max_step, bet, alf);
    
    % doing steps and updating active set:

    if is_aw
      % away step:
      alpha_t = (1+step)*alpha_t; % note that inactive should stay at 0;
      if abs(step - max_step) < 10*eps
          % drop step:
          number_drop = number_drop+1;
          alpha_t(id_A) = 0;
          I_active(I_active == id_A) = []; % remove from active set
          %TODO: could possibly also remove it from S_t
%           h = hashing(v_A);
%           mapping.delete(h);
%           
      else
          alpha_t(id_A) = alpha_t(id_A) - step;
      end
    else
      % FW step:
      alpha_t = (1-step)*alpha_t;

      % is this s_FW a new vertex?
      h = hashing(s_FW);
      if ~mapping.isKey(h)
          % we need to add vertex in S_t:
          max_index = max_index + 1;
          mapping(h) = max_index;
          S_t(:,max_index) = s_FW;
          id_FW = max_index;
          alpha_t(id_FW) = step; % this increase size of alpha_t btw
          I_active = [I_active, id_FW];
      else
          id_FW = mapping(h);
          if alpha_t(id_FW) < eps
              % we already had atom in 'correction poytope', but it was not
              % active, so now track it as active:
              I_active = [I_active, id_FW];
          end
          alpha_t(id_FW) = alpha_t(id_FW) + step;
      end

      % exceptional case: stepsize of 1, this collapses the active set!
      if step > 1-eps
          I_active = [id_FW];
      end
    end
    
    x_t = x_t + step * d; 
    
%     % remove zeros from active set if too many zeros in it
%     if sum(alpha_t == 0) > n
%         S_t = S_t(:, I_active);
%         alpha_t = alpha_t(I_active);
%         max_index = length(I_active);
%         I_active = 1:max_index;
%         % reconstruct mapping:
%         mapping = containers.Map();        
%         for index = 1:max_index
%             mapping(hashing(S_t(:,index))) = index; 
%         end
%     end
    
    % plot
    if doPlot
%         if mod(it, 1000) == 0
%             fprintf('iter %d/%d ', k, opts.N_it);
        if mod(it, 20) == 0
%             fprintf('.');
            [~,pp]=sort(x_t);
            imagesc(A(pp,pp)); title(sprintf('it %d',it));
            pause(0.1);
        end
    end

    [~, p_t] = sort(x_t);
    fPermVal = cost_fun(p_t);
    fPermVals(it-1) = fPermVal;
    if fPermVal < minf
        minf = fPermVal;
        bestPerm = p_t;
    end

    
end



res.primal = fvalues;
res.gap = gap_values;
res.number_away = number_away;
res.number_drop = number_drop;
res.S_t = S_t;
res.alpha_t = alpha_t;
res.x_t = x_t;
res.fPermVals = fPermVals;
res.bestperm = bestPerm;

% Return permutation
[~, p_t] = sort(x_t);


    % returns the id of the active atom with the worst value w.r.t. the
    % gradient
    function id = away_step(grad, S, I_active)
        s = grad' * S(:,I_active);
        [~,id] = max(s);
        id = I_active(id(1));
    end

     % MODIFY THE FOLLOWING FUNCTION for the atoms in your domain
    function h = hashing(permvec)
        % sequence is permutation vector
        % output is its conversion to a string
        h = char(join(string(permvec)));
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
        options.S_0 = x0;
        options.alpha_0 = 1;
        options.doPlot = false;
        options.A = eye(n);
        options.verbose = true;
    end

end