function [perm, x, res] = permManOpt(A, opts)


n = size(A,1);
% Gen. basis of the hyperplane that contains the permutation vectors
U = orthoConstBasis(n);
U = U(:,1:n-1);

opts_def = defaultOptions(n);
opts = build_opts(opts_def, opts);
manoptOpts = opts.manoptOptions;
if opts.dHuber == inf
    fh = @(x) two_SUM(x, A);
else
    % rescale dHuber by the radius of the sphere containing permutation
    % vectors
    xp = randperm(n)';
    yp = U'*xp;
    dilat = norm(yp);
    dHuber = opts.dHuber / dilat;
    fh = @(x) huberSUM(x, A, dHuber);
end

% Define function on the hyperplane
fhHpp = @(y)transfo2permhyperplane(U, y, fh);


% Define problem for manopt (Boumal etc.)
problem = struct;
manifold = spherefactory(n-1);
problem.M = manifold;
problem.costgrad = fhHpp;
y_0 = U'*opts.x_0;

% Numerically check gradient consistency (optional).
checkgradient(problem);

% Solve manifold problem
[y, ycost, info, manoptOpts] = trustregions(problem, y_0, manoptOpts);
x = U*y;
[~,perm] = sort(x);
res = struct;
res.ycost = ycost;
res.info = info;
res.manoptOpts = manoptOpts;


    function options = defaultOptions(n)
        options.x_0 = randperm(n)';
        options.manoptOptions = [];
        options.dHuber = inf;
    end
        

end


