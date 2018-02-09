
%%
p = 20;
sparam = p/2;
b = 1 + rand(p,1);
b = sort(b, 'descend');
bsum = sum(b);
a = 1.2*bsum;
elnet=0.3;
% Run convex programming software to perform projection
cvx_begin
    variable y(p);
    minimize elnet*norm(y - b, 1) + (1-elnet)*norm(y - b,2)
    subject to
        sum(y) == a;
        y >= 0;
cvx_end

myy = projconssparse(b,a,floor(p/2));

slop = [x,info] = Adlas(A,b,lambda,options)

figure; plot(y,'r+'); hold on; plot(myy, 'kx'); hold on; plot(b,'o'); %hold on; plot(myy1,'gd');