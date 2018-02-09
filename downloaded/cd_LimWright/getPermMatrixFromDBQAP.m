function [P,best_gap]=getPermMatrixFromDBQAP(D,Ap,S,nbTrial)
% get permutation matrix from doubly stochastic matrix (best of 'nbTrial' random trials)
n          = size(D,1);
% PermMatrix = eye(n);
% P          = PermMatrix;
u= (1:n)'+randn(n,1)/2;
v = D*u;
[~,perm] = sort(v);
[~,perm] = sort(perm);

PermMatrix = zeros(n);
PermMatrix(sub2ind([n,n],1:n,perm')) = 1;
P = PermMatrix;
% gaps       = zeros(nbTrial,1);

obj = get_cost(PermMatrix);

% nbTrial=size(y,2);

for i=1:nbTrial
    
%     u = cumsum( rand(n,1) );
%     u=y(:,i);
    
    u= (1:n)'+randn(n,1)/2;
    v = D*u;
    [~,perm] = sort(v);
    [~,perm] = sort(perm);

    PermMatrix = zeros(n);
    PermMatrix(sub2ind([n,n],1:n,perm')) = 1;

    newObj=get_cost(PermMatrix);
    if  newObj < obj
        obj=newObj;
        P=PermMatrix;
    end  
%     gaps(i) = obj - get_cost(D);
end




%[max(gaps),min(gaps),mean(gaps),std(gaps)]

if ( [1,zeros(1,n-1)]*P*(1:n)'  >=  [zeros(1,n-1), 1]*P*(1:n)' )
	P = P(:,n:-1:1);
end


best_gap = obj - get_cost(D);

	% nested function
	function cost=get_cost(arg)
        cost = sum(dot(-S,arg'*Ap*arg));
		%cost = trace(y'*arg'*L*arg*y) - mu*norm(arg*y,'fro')^2;
% 		cost = trace(y'*arg'*L*arg*y) - mu*norm( (eye(n) - 1/n*ones(n))*arg,'fro')^2;
	end

end
