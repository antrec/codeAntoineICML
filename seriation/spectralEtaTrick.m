function [perm, fs] = spectralEtaTrick(A, opts)
% Uses the variational trick |x| = argmin_{y} 1/y * x^2 + y to solve the
% 1SUM problem [sum_ij A_ij |x_i - x_j|, x \in \Permutations] or Huber
% version [sum_ij A_ij hub(x_i - x_j, delta), with hub(x,delta) = min(x^2,
% delta*(2*|x| - delta))]

n = size(A,1);

% Build options
opts_def = defaultOptions();
opts = build_opts(opts_def, opts);
dh = opts.dHuber;
fh = @(x) huberSUM(x,A,dh);
doAdd = false;
if opts.addIter
    doAdd = true;
    gam = opts.addGamma;
end
useFid = opts.useFiedler;
doPlot = opts.doPlot;
maxNoDecrease = opts.maxNoDecrease;

% Initialization
bestscore = +inf;
fs = zeros(opts.Niter,1);


% Initialize
if issparse(A)
    [is, js, vs] = find(A);
    etas = ones(size(vs));
else
    etas = ones(n);
    if ~useFid
        T = max(dh, abs(repmat((1:n)',1,n) - repmat((1:n),n,1)));
    end
end

noDecreaseCpt = 0;

for iter = 1:opts.Niter

    if issparse(A)

        % Approximately solve sum A_ij*(1/eta_ij*(x_i - x_j)^2 + eta_ij)
        myA = sparse(is, js, vs./etas, n, n);
        [mypi, myx] = spectralOneCC(myA);
        [~,myp] = sort(mypi);

        % Update eta and compute objective
        if useFid
            newetas = max(dh, abs(myx(is)-myx(js)));
            obj = fh(myx);
        else
            newetas = max(dh,abs(myp(is)-myp(js)));
            obj = fh(myp);
        end

    else

        % Approximately solve sum A_ij*(1/eta_ij*(x_i - x_j)^2 + eta_ij)
        myA = A./etas;
        [mypi, myx] = spectralOneCC(myA);
        [~,myp]=sort(mypi);

        % Update eta and compute objective
        if useFid
            newetas = max(dh, abs(repmat(myx',1,n) - repmat(myx,n,1)));
            obj = fh(myx);
        else
            newetas = T(myp,myp);
            obj = fh(myp);
        end

    end

    if doAdd
        etas = gam*newetas + (1-gam)*etas;
    else
        etas = newetas;
    end

    % Keep objective and update best permutation seen so far
    fs(iter) = obj;
    if obj < bestscore
        bestscore = obj;
        perm = mypi;
        noDecreaseCpt = 0;
    else
        noDecreaseCpt = noDecreaseCpt + 1;
    end
    
    if noDecreaseCpt > maxNoDecrease
        break;
    end
    
    % Plot
    if doPlot
        if issparse(A)
            subplot(1,2,1); spy(A(mypi,mypi)); colorbar;
%             subplot(1,3,2); scatshowsp(sparse(is,js,etas,n,n)); colorbar;
        else
            subplot(1,2,1); imagesc(A(mypi,mypi)); colorbar;
%             subplot(1,2,2); imagesc(etas); colorbar;
        end
        subplot(1,2,2); plot(fs(1:iter),'-d');
        title(sprintf('it%d, score (best) : %1.2e (%1.2e)',iter, obj, bestscore));
        pause(0.1);
    end

end
   

    function options = defaultOptions()
    % Default options
    options = struct;
    options.doRegLap = false;
    options.dHuber = 1;
    options.useFiedler = false;
    options.Niter = 30;
    options.addIter = true;
    options.addGamma = 0.5;
    options.doPlot = false;
    options.maxNoDecrease = 10;

    end

end