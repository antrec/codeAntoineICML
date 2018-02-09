function [x, active_set, alphas, scores] = permFW(func_h, n, opts, S)
% Frank-Wolfe in Permutathedron

% parse options and sets undefined ones to default
opts_def.N_it = 1e3; % Max number of iterations
opts_def.eps = 1e-3;
opts_def.FW_type = 'reg'; % regular (reg) or away-step (away) FW
opts_def.ls_alpha = 1e-4;
opts_def.ls_beta = 0.9;
opts_def.x_i = randperm(n)';
opts_def.w_bias = zeros(n,1);
opts_def.keep_active = false;
if nargin <= 2
    opts = opts_def;
else
    opts = build_opts(opts_def, opts);
end

alpha = opts.ls_alpha;
beta = opts.ls_beta;
scores = zeros(1, opts.N_it);
active_set = zeros(n, opts.N_it);
active_set(:,1) = opts.x_i;
alphas = zeros(1, opts.N_it);
alphas(1)=1;
% active_set = [opts.x_i];
% alphas = [1];
lprec = -1e10;
lb = zeros(1,opts.N_it);
c = 1./2*(n+1)*ones(n,1);
x = 0.5*(c + opts.x_i);
w_bias = opts.w_bias;

for k=1:opts.N_it
    
    % compute objective and gradient
    [obj, grad] = func_h(x);
%     ff = func_h(x);
%     obj = ff.obj;
    scores(k) = obj;
%     grad = ff.grad;

    % compute descent direction
%     s = LMO_sort_with_tiebreak_modif(grad);
    [~,ss]=sort(-grad);
    [~,s] = sort(ss);
%     s = proj_half_phthdrn(s, w_bias);
    dfw = s - x; % FW descent direction
    
    if strcmp(opts.FW_type, 'away')
        
        % compute away (ascent) direction
        asidx = alphas>0;
        as = active_set(:, asidx);
        alph = alphas(asidx);
        activeLinScore = as'*grad;
        [~,ai] = max(activeLinScore);
        v = as(:,ai);
        da = x - v; % away direction

        %chose between FW and away step according to gap
        if -grad'*dfw >= -grad'*da % FW step
            dt = dfw;
            gamMax = 1;
            stepType = 'fw';
        else % away step
            dt = da;
            gamMax = alph(ai)./(1-alph(ai));
            stepType = 'away';
        end

        % line search for generic convex objective
        gamm = gamMax;
%         ff = func_h(x);
%         fx = ff.obj;
        fx = func_h(x);
        fxx = func_h(x + gamm*dt);
%         ff = func_h(x + gamma*dt);
%         fxx = ff.obj;
        
        while (fxx > fx + alpha*gamm*grad'*dt)
            gamm = beta*gamm;
            fxx = func_h(x + gamm*dt);
%             ff = func_h(x + gamma*dt);
%             fxx = ff.obj;
        end
        
        if strcmp(stepType,'fw') % Frank-Wolfe step
            if gamm == 1
                alphas(1:k-1) = 0;
                active_set(:,k) = s;
                alphas(k) = 1;
%                 active_set = s;
%                 alphas = 1;
            elseif gamm >0
                [Lia, Locb] = ismember(s',as','rows');
                if Lia
                    alph = (1-gamm)*alph;
                    alph(Locb) = alph(Locb) + gamm;
                    alphas(asidx) = alph;
                else
                    active_set(:,k) = s;
                    alphas(k) = gamm;
                    alphas(asidx) = (1-gamm)*alph;
%                     active_set = [active_set, s];
%                     alphas = (1-gamma)*alphas;
%                     alphas = [alphas, gamma];
                end
            end

        else % away-step
            if gamm == gamMax % drop step
                fasidx = find(asidx);
                aai = fasidx(ai);
                alphas(aai) = 0;
%                 active_set(:,aai) = zeros(n,1);
%                 active_set = [active_set(:,1:ai-1), active_set(:,ai+1:end)];
%                 alphas = [alphas(1:ai-1), alphas(ai+1:end)];
                alphas = (1+gamm).*alphas;
            else
                alphas = (1+gamm)*alphas;
                fasidx = find(asidx);
                aai = fasidx(ai);
                alphas(aai) = alphas(aai) - gamm;
            end
        end    
        
    else % "simple" FW with no away steps
        
        gamm = 2./(2+k); stepType='fw'; dt = dfw;
        
        % if one wants to keep active set (slower)
        if opts.keep_active
            
            asidx = alphas>0;
            as = active_set(:, asidx);
            alph = alphas(asidx);
            [Lia, Locb] = ismember(s',as','rows');
            if Lia
                alph = (1-gamm)*alph;
                alph(Locb) = alph(Locb) + gamm;
                alphas(asidx) = alph;
            else
                active_set(:,k) = s;
                alphas(k) = gamm;
                alphas(asidx) = (1-gamm)*alph;
%                     active_set = [active_set, s];
%                     alphas = (1-gamma)*alphas;
%                     alphas = [alphas, gamma];
            end
            
        end
        
    end
    
    % update x
    x = x + gamm*dt;

    % keep lower bound
    lb(k) = max(lprec, obj + dt'*grad);
    lprec = lb(k);
    
    % print iteration
    if mod(k, 1000) == 0
        fprintf('iter %d/%d ', k, opts.N_it);
    elseif mod(k, 100) == 0
        fprintf('.');
        [~,pp]=sort(x);
        imagesc(S(pp,pp)); title(sprintf('it %d',k));
        pause(0.001);
    end
    
    % remove zeros from active set
    active_set = active_set(:, alphas>0);
    alphas = alphas(alphas>0);
        
end
