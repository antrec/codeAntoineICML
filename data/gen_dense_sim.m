function [simMatrix] = gen_dense_sim(n, noise, alpha)
% simulate pairwise comparisons based on true ordering (=1:n), i.e.
% choose randomly inequality of type (PI_true*omega)(i)<(PI_true*omega)(j) to
% incorporate as additional info
if nargin <= 2
    alpha=1;
end
if nargin == 1
    noise=0.;
end
% p=alpha*n*(n-1)/2;
% p=10*n*log(n);
% p=50*n; %number of (different) pairwise comparisons
omega=(1:n)';
permOmega=omega;
% permOmega=randperm(n)';
truePI = zeros(n);
truePI(sub2ind([n,n], permOmega', 1:n)) = 1;
pairwiseMat=triu(ones(n)) - tril(ones(n),-1) ;
ratioNbComp=alpha;
triIdx=find(triu(ones(n),1));
isObserved=zeros(n); 
isObserved(triIdx)=rand(n*(n-1)/2,1)<=ratioNbComp;
isObserved=isObserved + isObserved';
pairwiseMat(~isObserved)=0;
isNoisy=zeros(n); 
isNoisy(triIdx)=rand(n*(n-1)/2,1)<=noise;
isNoisy=isNoisy + isNoisy';
pairwiseMat(logical(isNoisy))=-pairwiseMat(logical(isNoisy));
pairwiseMat(logical(eye(n)))=1;

pairwiseMat=truePI'*pairwiseMat*truePI;
% isObserved=truePI'*isObserved*truePI;

% construct similarity matrix based on pairwise comparisons
simMatrix=pairwiseMat*pairwiseMat' + n ;
% nbMatch=isObserved;
% simMatrix=zeros(n);
% for i=1:n
%     for j=1:n
%         for k=1:n
%             if nbMatch(i,k)>0 && nbMatch(j,k)>0
%                 simMatrix(i,j)=simMatrix(i,j) +1-1/2*abs(pairwiseMat(i,k)-pairwiseMat(j,k));
%             else
%                 simMatrix(i,j)=simMatrix(i,j)+1/2;
%             end
%         end
%     end
% end
simMatrix=simMatrix-min(simMatrix(:));