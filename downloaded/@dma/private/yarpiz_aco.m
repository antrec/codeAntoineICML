%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YOEA103
% Project Title: Ant Colony Optimization for Quadratic Assignment Problem
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

%Slightly modified and repackaged to make it callable as a function, 2016.
%It minimises the QAP problem <P*D*P',W> where D and W are assumed symmetric.
function [perm, details] = yarpiz_aco(D, W, MaxIt)

%D = max(max(D))-D;
%W = -W;

model.w = W;
model.d = D;
model.n = size(D,1);
model.m = size(D,1);

CostFunction=@(p) MyCost(p,model); 
nVar=size(D,1);     %Number of Decision Variables



%ACO Parameters
%MaxIt=500;      % Maximum Number of Iterations

nAnt=50;        % Number of Ants (Population Size)
Q=1;
tau0=10;        % Initial Phromone
alpha=0.3;      % Phromone Exponential Weight
rho=0.1;        % Evaporation Rate


% Initialization
tau=tau0*ones(model.m,nVar);

BestCost=zeros(MaxIt,1);    % Array to Hold Best Cost Values

% Empty Ant
empty_ant.Tour=[];
empty_ant.Cost=[];

% Ant Colony Matrix
ant=repmat(empty_ant,nAnt,1);

% Best Ant
BestSol.Cost=inf;


% ACO Main Loop

for it=1:MaxIt
    
    % Move Ants
    for k=1:nAnt
        
        ant(k).Tour=[];
        
        for l=1:nVar
            
            P=tau(:,l).^alpha;
            
            P(ant(k).Tour)=0;
            
            P=P/sum(P);
            
            j=RouletteWheelSelection(P);
            
            ant(k).Tour=[ant(k).Tour j];
            
        end
        
        ant(k).Cost=CostFunction(ant(k).Tour);
        
        if ant(k).Cost<BestSol.Cost
            BestSol=ant(k);
        end
        
    end
    
    % Update Phromones
    for k=1:nAnt
        
        tour=ant(k).Tour;
        
        for l=1:nVar
            
            tau(tour(l),l)=tau(tour(l),l)+Q/ant(k).Cost;
            
        end
        
    end
    
    % Evaporation
    tau=(1-rho)*tau;
    
    % Store Best Cost
    BestCost(it)=BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
   
    perm = (BestSol.Tour);
    details = BestCost;
    
end

end


%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YOEA103
% Project Title: Ant Colony Optimization for Quadratic Assignment Problem
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
function z=MyCost(p,model)

    %n=model.n;
    w=model.w;
    d=model.d;

   % z=0;
    %for i=1:n-1
    %    for j=i+1:n
    %       z=z+w(i,j)*d(p(i),p(j));
    %        
     %   end
    %end
    
    z = sum(sum(w.*d(p,p))); %modified to be faster

end



%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YOEA103
% Project Title: Ant Colony Optimization for Quadratic Assignment Problem
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function j=RouletteWheelSelection(P)

    r=rand;
    
    C=cumsum(P);
    
    j=find(r<=C,1,'first');

end



