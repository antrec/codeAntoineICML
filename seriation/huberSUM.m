function [f, df] = huberSUM(x, A, dh)
% Compute objective and gradient of the huberized version of the 2SUM
% objective, namely sum_{ij} A_{ij} hub(x_i - x_j, dh), where hub(x,dh) = 
% min(x^2,delta*(2*|x| - delta))

    n = size(A, 1);
    
    doTrans=false;
    if size(x) ~= [n,1]
        if size(x) == [1,n]
            x = x';
            doTrans = true;
        else
            fprintf('x must be of length n !');
        end
    end
    
    if issparse(A)

        [ii, jj, vv] = find(A);

        dx = x(ii) - x(jj);
        mask = (abs(dx) > dh);

        dxh = zeros(size(dx));
        dxh(~mask) = dx(~mask).^2;
        dxh(mask) = 2*dh*(abs(dx(mask)) - 1/2*dh);

        dg = zeros(size(dx));
        dg(~mask) = 2*dx(~mask);
        dg(mask) = 2*dh*sign(dx(mask));

        objA = sparse(ii, jj, vv.*dxh, n, n);
        gradA = sparse(ii, jj, vv.*dg, n, n);

    else

        dX = repmat(x, 1, n) - repmat(x', n, 1);
        mask = abs(dX)>dh;

        dXh = zeros(size(dX));
        dXh(~mask) = dX(~mask).^2;
        dXh(mask) = 2*dh*(abs(dX(mask)) - 1/2*dh);

        dG = zeros(size(dX));
        dG(~mask) = 2*dX(~mask);
        dG(mask) = 2*dh*sign(dX(mask));

        objA = A.*dXh;
        gradA = A.*dG;

    end

    f = sum(sum(objA));
    df = sum(gradA, 2);

    if issparse(A)
        f = full(f);
        df = full(df);
    end
    
    if doTrans
        df = df';
    end

end
