
%%
p = 20;
sparam = p;%/2;
b = 1 + rand(p,1);
b = sort(b, 'descend');
b(end-5:end)=0;
bsum = sum(b);
a = 1.5*bsum;
elnet=0.3;
% Run convex programming software to perform projection
cvx_begin
    variable y(p);
    minimize elnet*norm(y - b, 1) + (1-elnet)*norm(y - b,2)
    subject to
        sum(y) == a;
        y >= 0;
cvx_end

sparam = min(sparam, nnz(b));
myy = projconssparse(b,a,sparam);
figure; plot(y,'r+'); hold on; plot(myy, 'kx'); hold on; plot(b,'o'); %hold on; plot(xslope,'gd');
%%
lambdas = 1./(1:p).^2;
sumslope=-inf;
mu = 1e-2;
while sumslope < a*(1-1e-3)
    mu = 1.3*mu;
%     lambdas = 0.5*lambdas;
% lambdas = ones(
    [xslope,info] = Adlas(eye(p),b,lambdas);
    sumslope = sum(xslope);
end
%%

figure; plot(y,'r+'); hold on; plot(myy, 'kx'); hold on; plot(b,'o'); %hold on; plot(xslope,'gd');