function [perm] = seriationRobGM(A, opts)

n = size(A,1);
opts_def = defaultOptions(n);
opts = build_opts(opts_def, opts);
T = abs(repmat((1:n)',1,n) - repmat((1:n),n,1));
iter = opts.iter;
p0 = opts.p_0;
paramgm = struct;
paramgm.maxIter = iter;
paramgm.verbose = 1;
paramgm.P0 = 0.9*eye(n) + 0.1*ones(n);

switch opts.Toeplitz
    case '2SUM'
        T = T.^2;
    case 'R2SUM'
        dh = opts.dH;
        T = min(dh^2, T.^2);
    case 'Huber'
        dh = opts.dH;
        T = min(T.^2, dh*(2*T - dh));
    otherwise
        fprintf('invalid opts.Toeplitz. Solving 1SUM');
end

% Scale T with A
nA = norm(A,'fro');
nT = norm(T,'fro');
T = nA/nT * T;

% Check if just QAP for Robust Seriation works
[~,Pp]=graph_matching(max(T(:))-T,A(p0,p0), paramgm);
% [~,Pp]=graph_matching(max(T(:))-T,A, paramgm);

perm = Pp*(1:n)';
perm = p0(perm);
% res.f = f;
% res.x = x;
% res.fs = fs;
% res.myps = myps;



function options = defaultOptions(n)
    options.p_0 = randperm(n)';
    options.Toeplitz = '2SUM';
    options.dH = floor(n/10);
    options.iter = 10000;
end

end