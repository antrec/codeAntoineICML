%From paper: A Graduated Non-Convexity Relaxation for Large Scale Seriation.
%            X. Evangelopoulos, A.J. Brockmeier, T. Mu and J.Y. Goulermas, SIAM Data Mining, 2017.
%
%D : an nxn symmetric dissimilarity matrix
%mu_init : initial regularization parameter (\mu_1 in the paper)
%gamma : continuation rate (>1)
%gnc_perm : final permutation
%details.two_sum : final 2-SUM score
%
%X. Evangelopoulos, 2016.

function [gnc_perm, details] = gncr_twosum(D,varargin)

A = max(max(D)) - D; %convert to similarity for this method
n = size(A,1);

%Laplacian
lapfun = @(x) (diag(sum(x,2)) - x);
L_A = lapfun(A);

%lambda values to control convex and concave regularization
lambda2 = eigs(L_A,2,'sa');
lambda2 = max(lambda2);
lambdaN = eigs(L_A,1,'la');


%default parameters
mu_init = lambda2;
gamma = 1.05;

if mod(length(varargin),2)
    error('dma:gncr', 'missing parameters')
end

for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
        case 'mu_init'
            mu_init = varargin{i+1};
        case 'gamma'
            gamma = varargin{i+1};
        case 'delta'
        otherwise
            error('dma:gncr', 'unknown parameter: %s', varargin{i} )
    end
end

e = (1:n)';
nstep = 500; %number of FW steps
H = eye(n) - (1/n)*ones(n); %centering matrix

%Conditional Gradient approach
perm = e;
gperm = (e+(n+1)/2*ones(n,1))/2; %starting point
mu = mu_init;

while mu <= gamma*lambdaN
    for j = 1:nstep
        
        %compute conditional gradient & instead of lap do sort
        [~,perm] = sort(L_A*gperm-mu*H*gperm,'descend');
        [~,iperm] = sort(perm,1);
        gperm_1 = iperm;
        
        %Line search
        c1 = gperm_1'*(L_A-mu*H)*gperm_1;
        c2 = gperm'*(L_A-mu*H)*gperm;
        c3 = gperm'*(L_A-mu*H)*gperm_1;
        
        if c1+c2-2*c3 <= 0
            alpha = 1;
            
            if c1 >= c2
                break;
            end
        else
            alpha = (c2-c3)./(c1+c2-2*c3);
            alpha = max(0,min(1,alpha));
            
            %Convergence check
            if alpha <= 1e-2
                break;
            end
        end
        
        gperm = alpha*gperm_1+(1-alpha)*gperm;
        
        %get permutation of next step
        [~,perm] = sort(gperm);
        
    end
    
    %mu update
    if mu == 0 %handle the case where mu_init == 0
        mu = 0.001;
    end
    mu = mu*gamma;
    
end
W2 = abs(dma.template(n,'squared'));

gnc_perm = perm;
details.two_sum = trace(D(gnc_perm,gnc_perm)*W2);
end


