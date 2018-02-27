%From paper: Approximation Methods for Large Scale Object Sequencing.
%            X. Evangelopoulos, A.J. Brockmeier, T. Mu and J.Y. Goulermas,
%            Machine Learning, 2018. (under submission)
%
%D : an nxn symmetric dissimilarity matrix
%delta :  pseudo-Huber parameter
%mu_init : initial regularization parameter (\mu_1 in the paper)
%gamma : continuation rate (>1)
%hgnc_perm : final permutation
%details.one_sum : final 1-SUM score
%
%X. Evangelopoulos, 2017.


function [hgnc_perm, details] = gncr_pseudo_huber(D, varargin)

if mod(length(varargin),2)
    error('dma:gncr', 'missing parameters')
end

for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
        case 'delta'
            delta = varargin{i+1};
    end
end

 if ~exist('delta', 'var')
   error('dma:gncr', 'missing parameter delta')
 end

 
 
A_rp = max(max(D)) - D; %convert to similarity for this method
n = size(A_rp,1);

e = (1:n)';
nstep = 500; %number of FW steps
H = eye(n) - (1/n)*ones(n); %centering matrix

gperm = (e+(n+1)/2*ones(n,1))/2; %starting point

%Hessian \phi_delta(x)
[~,~,Hess] = huberhess(gperm,A_rp,delta);

%lambda values for convex and concave regularization
lambdas = sort(eig(Hess));
lambda2 = lambdas(2);

%default parameters
mu_init = lambda2;
gamma = 1.05;

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


%Conditional Gradient approach
mu = mu_init;

while ~isequal(sort(gperm),e)% increase mu until we get a discrete solution
    for j = 1:nstep
        W=bsxfun(@minus,gperm,gperm');
        
        %compute conditional gradient & instead of lap do sort
        grad = (1/delta)*sum(A_rp.*(W./(sqrt(1+(W.^2)./delta^2))),2);
        [~,perm]=sort(grad-mu*H*gperm,'descend');
        [~,iperm]=sort(perm,1);
        gperm_1 = iperm;
        
        %Line search
        gperm_a =@(alpha) (gperm_1*alpha +gperm*(1-alpha));
        huber_w =@(alpha) -2*(gperm_a(alpha)*gperm_a(alpha)')...
            + bsxfun(@plus, gperm_a(alpha).^2,gperm_a(alpha).^2');
        eval = @(alpha) A_rp(:)'*(delta*(sqrt(1+(reshape(huber_w(alpha),[],1)/delta^2)))-delta) -mu*gperm_a(alpha)'*H*gperm_a(alpha);
        
        options = optimset('Display','off','MaxFunEvals',500,'MaxIter',500);
        [alpha,~,~,~]=fminbnd(eval,0,1,options);
        
        gperm = gperm_a(alpha);
        
        %Convergence check
        if alpha < 1e-2 || eval(alpha)==eval(0)
            if eval(0) == eval(1)
               gperm = gperm_1;
            end
            break;
        end        
    end
     %mu update
     if mu == 0 %handle the case where mu_init == 0
         mu = 0.001;
     end
    mu = mu*gamma;
end
W1 = abs(dma.template(n,'linear'));%absolute positional differences template

[~,hgnc_perm] = sort(gperm);
details.one_sum = trace(D(hgnc_perm,hgnc_perm)*W1);
end



%--------------------------------------------------------------------------------------------------




%Helper function to pshuber_gncr.m - Computes the Hessian \phi_delta(x) at 
%x, for a given similarity matrix and a given delta, as defined in Eq. (22)
%in the paper (Evangelopoulos et al., 2018).
%
%A_rp : input similarity matrix
%delta: input delta parameter
%fun : objective function value at x
%grad : gradient at x
%H : Hessian at x
%
%X. Evangelopoulos, 2017.
function [fun,grad,H] = huberhess(x,A_rp,delta)

W = bsxfun(@minus,x,x'); 
        
g_f = (1/delta)*sum(A_rp.*(W./(sqrt(1+(W.^2)./delta^2))),2);

huber_w = -2*(x*x') + bsxfun(@plus, x.^2,x.^2');

grad = 2*g_f;

Hess = (2/delta)*(A_rp./(sqrt(1+(W.^2)./delta^2).^3));

n = length(x);

Hess(1:n+1:n*n) = 0; %offdiagonal elements of the Laplacian (G)

Hess_diag = sum(Hess,2); %diagonal elements of the Laplacian (diag(G1))

H = diag(Hess_diag) - Hess; %L_G = diag(G1) - G

fun = A_rp(:)'*(delta*(sqrt(1+(reshape(huber_w,[],1)/delta^2)))-delta);
end