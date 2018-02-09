function [perm, res] = seriationLWCD(A, opts)

n = size(A,1);
opts_def = defaultOptions(n);
opts = build_opts(opts_def, opts);
T = abs(repmat((1:n)',1,n) - repmat((1:n),n,1));
iter = opts.iter;
p0 = opts.p_0;

switch opts.Toeplitz
    case '2SUM'
        T = T.^2;
    case 'R2SUM'
        dh = opts.dH;
        T = min(dh^2, T.^2);
    case 'Huber'
        dh = opts.dH;
        T = min(T.^2, dh*(2*T - dh));
    % otherwise T is unchanged --> 1SUM
end
        
% Check if just QAP for Robust Seriation works
[f, perm] = cd_process(T, A(p0,p0), iter, opts.heur, opts.reg_limit);
perm = p0(perm);
res.f = f;



function options = defaultOptions(n)
    options.p_0 = randperm(n)';
    options.Toeplitz = '2SUM';
    options.dH = floor(n/10);
    options.iter = 300;
    options.heur = 'full';
    options.reg_limit = 30;
end

end