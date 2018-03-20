function [Zt, St] = seriationDuplialtProj(A, c, RobSerFunc, opts, ubval)
% alternate projection to solve seriation with duplications


% t = clock;
nsmall = size(A,1);
[n, opts_def] = defaultOpts(A,c);
opts = build_opts(opts_def, opts);

Ztrue = opts.Ztrue;
Zt = opts.Zi;
dcm = diag(1./c);
St = Zt'*dcm*A*dcm'*Zt;
doPlot = opts.doPlot;
Niter = opts.Niter;
doQuant = opts.doQuantileNorm;
k = opts.ksparse;
dh = opts.dh;
doWeakR = opts.doWeakR;

optsProj = opts.optProj;
optsProj.dh = dh;
optsProj.ub = max(ubval);

% test something (quantile normalization?)
if doQuant
    onestri = tril(ones(n),0);
    onestri = onestri(:);
    nnzitri = find(onestri);
    % Use only lower triangle to keep matrices symmetric
    svals = sort(opts.Svalues(nnzitri));
end

% et = etime(clock, t);
% fprintf('read options in %3.3fs ', et);
% t = clock;
% Build projection on R matrix constraints here
[square2ltmat, lineqmat] = buildRConsSparse(n, k);
% et = etime(clock, t);
% fprintf('computed R-projection matrices in %3.3fs ', et);


iter = 100;
x0 = eye(n);
T = max(0,abs(repmat((1:n)',1,n) - repmat((1:n),n,1)));


for it=1:Niter
    
    perm = spectralOneCC(St);
    if it > 1
        % Initialize with spectral
%         pspectr = spectralOneCC(St);
         pspectr = perm;
        % Run Robust Seriation
        perm = RobSerFunc(St(pspectr,pspectr));
        perm = pspectr(perm);
%         perm = pspectr(perm);

        if size(perm,1) < size(perm,2)
            perm = perm';
        end
    end
    
    % Perform projection on pre-R
    t = clock;
%     [Sproj] = proj2RmatSparse(St(perm,perm), k, square2ltmat, lineqmat);
%     if mod(it,5) == 0
%         doWeakR = true;
%     else
%         doWeakR = false;
%     end

    if doPlot
        subplot(2,2,4); imagesc(St(perm,perm));
%         plot(perm,'o');
    end

%     if it > Niter/2 && opts.doWeakR
%         doWeakR = true;
%     else
%         doWeakR = false;
%     end
%     [Sproj] = proj2RmatWithLinProg(St(perm,perm),doWeakR, ubval);
    projRopts = [];
    projRopts.ubs = ubval;
    Sproj = proj2Rmat2(St(perm,perm),projRopts);
%     Sproj = proj2RmatAll(St(perm,perm));
    et = etime(clock, t);
    n_put = sum(perm == (1:n)');
%     fprintf('R-proj in %1.2fs, %d unchanged elements.\n', et, n_put);
%     pause(1);
    if n_put == n
        break;
    end
    
    Zt = Zt(:, perm);

    % test something (quantile normalization?)
    if doQuant
%         Sp2 = gen_dense_sim(n,0,1);
%         Sproj = sum(sum(Sproj))/sum(sum(Sp2)) * Sp2;
        myvals = Sproj(:);
        myvals = myvals(nnzitri);
        [~,idxs] = sort(myvals);
        [~,idxs]=sort(idxs);
        mynormval = zeros(n^2,1);
        mynormval(nnzitri) = svals(idxs);
%         mynormval = svals(idxs);
%         zerosidx = find(myvals==0);
%         mynormval(nnzitri(zerosidx))=0;

        Sproj = reshape(mynormval, n, n); % ! A.R. this is a triangular matrix
    end
    t = clock;

    
%     Sproj = tril(Sproj) - tril(Sproj,-dh);
    St = proj2DuplConsSparseClosedForm(Sproj, Zt, A, optsProj);

    if doPlot
        subplot(2,2,1); imagesc(St); colorbar; title('S(t)');
        
            transpDists = zeros(2,nsmall);
            Zinv = Zt(:,end:-1:1);
            for ii=1:nsmall
                transpDists(1,ii) = (1/sum(Ztrue(ii,:)))* eval_twins(full(Ztrue(ii,:)),full(Zt(ii,:)));
                transpDists(2,ii) = (1/sum(Ztrue(ii,:)))* eval_twins(full(Ztrue(ii,:)),full(Zinv(ii,:)));
            end
            meanDists = mean(transpDists(1,:));
            meanInv = mean(transpDists(2,:));
            if meanInv < meanDists
                subplot(2,2,2); spy(Zinv,'.'); hold on; spy(Ztrue, 'rd'); hold off;
                title(sprintf('it %d | %1.2f', it, meanInv));
            else
                subplot(2,2,2); spy(Zt,'.'); hold on; spy(Ztrue, 'rd'); hold off;
                title(sprintf('it %d | %1.2f', it, meanDists));
            end
%         subplot(2,2,4);% imagesc(St(perm,perm));
%         plot(perm,'o');
        subplot(2,2,3); imagesc(Sproj); colorbar;
        pause(0.1);
%         pause;
    end
    
    % reduce number of constraints for R projection
%     [ii,jj,~]=find(St);
%     k = max(abs(ii-jj));
%     [square2ltmat, lineqmat] = buildRConsSparse(n, k);

end

    function [n,options] = defaultOpts(A,c)
        smalln = length(c);
        n = sum(c);
        iis = []; jjs = []; cpt=0;
        for i=1:smalln
            iis = [iis, i*ones(1, c(i))];
            jjs = [jjs, cpt+(1:c(i))];
            cpt = cpt + c(i);
        end
        options.Zi = sparse(iis, jjs, ones(1, length(iis)), smalln, n);
        options.Ztrue=options.Zi;
        options.doPlot = false;
        options.Niter = 50;
        options.doQuantileNorm = false;
        dc = diag(1./c);
        Si = (options.Zi)'*dc*A*dc'*(options.Zi);
        options.Svalues = Si(:);
        options.ksparse = n;
        options.dh = n;
        optionsProj = struct;
        optionsProj.dh = options.dh;
        optionsProj.softCons = false;
        optionsProj.sparsprop = 1;
        options.optProj = optionsProj;
        options.doWeakR = false;
    end


end