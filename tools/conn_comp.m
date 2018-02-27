function [conComp,nbComp,singles]=conn_comp(S)
[labels, nbComp] = graph_conn_comp(S);%    conncomp(S);
% [nbComp, labels] = graphconncomp(S);
% [labels, nbComp] = conncomp(S);

conComp=cell(nbComp,1);
lengthComp=zeros(1,nbComp);
for i=1:nbComp
    conComp{i}=find(labels==i);
    lengthComp(i)=numel(conComp{i});
end
% reorder contigs according to their size
[~, idx]=sort(lengthComp,'descend');
conComp=conComp(idx);
lengthComp=lengthComp(idx);
firstSingle=find(lengthComp==1,1,'first');
singles=zeros(1,nbComp-firstSingle+1);
for i=firstSingle:nbComp
    singles(i-firstSingle+1)=conComp{i};
end