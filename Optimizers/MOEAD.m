% function [Population,FunctionValue,Archive]=MOEAD(Problem,Generations,num_dimensions,Population,N)
    
function [Population,FunctionValue,Archive]=MOEAD(Problem,MaxFun,num_dimensions,Population,N)

%     Solutions_gen = struct;
    addpath(genpath('MOEAD_files'));
    addpath(genpath('Public'));
    M = Problem.numberOfObjectives;
%     if  M ==2
%         N = 100; % population size
%     elseif M ==3
%         N = 150;
%     else
%         N = 300;
%     end

    MaxValue = 1*ones(1,num_dimensions);
    MinValue = -1*ones(1,num_dimensions);
    Boundary = [MaxValue;MinValue];
%     D = size(Boundary,2);
%     Population = rand(N,D);
%     Population = repmat(MinValue,N,1) + (repmat(MaxValue,N,1) - repmat(MinValue,N,1)).*rand(N,num_dimensions);
%     Population = Population.*repmat(MaxValue,N,1)+(1-Population).*repmat(MinValue,N,1); % generating random population
    
    FunctionValue = zeros(size(Population,1),M);
    for i = 1:size(Population,1)
        FunctionValue(i,:) = Problem.evaluate(Population(i,:)); 
    end
    Archive = [Population,FunctionValue]; % storing all evaluated solutions
%     Solutions_gen(1).data = [Population,FunctionValue]; 
    W = ref_vectors(M);
    N = size(W,1);
    T = ceil(N/1);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %% Generate random population
%     Population = Global.Initialization();
    Z = min(FunctionValue,[],1);

    nfun = 0;
    %% Optimization
    while nfun <= MaxFun   
%         Gene
        % For each solution
        for i = 1 : N    
            % Choose the parents
            P = B(i,randperm(size(B,2)));

            % Generate an offspring
%             MatingPool = F_mating(Population,N); % mating pool
            Coding = 'Real';
%             Offspring = P_generator(MatingPool,Boundary,Coding,N);
            MatingPool = Population(P(1:2),:);
            Offspring = P_generator(MatingPool,Boundary,Coding,1);
            Fitness = zeros(size(Offspring,1),M);
            
            nfun = nfun + size(Offspring,1);
    
            for o = 1:size(Offspring,1)
                Fitness(o,:) = Problem.evaluate(Offspring(o,:));
            end
%             Archive = [Archive;[Offspring,Fitness]];
            % Update the ideal point
            Z = min(Z,Fitness);

            % Update the neighbours
            
             
                    % PBI approach
%             normW   = sqrt(sum(W(P,:).^2,2));
%             normP   = sqrt(sum((FunctionValue(P,:)-repmat(Z,T,1)).^2,2));
%             normO   = sqrt(sum((Fitness-Z).^2,2));
%             CosineP = sum((FunctionValue(P,:)-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
%             CosineO = sum(repmat(Fitness-Z,T,1).*W(P,:),2)./normW./normO;
%             g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
%             g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
            
            g_old = max(abs(FunctionValue(P,:)-repmat(Z,T,1)).*W(P,:),[],2);
             g_new = max(repmat(abs(Fitness-Z),T,1).*W(P,:),[],2);
            
            r = g_new<=g_old;
            t = P(r);
            Population(t,:) = repmat(Offspring,length(t),1);
            FunctionValue(t,:) = repmat(Fitness,length(t),1);
            
%             Population(P(r),:) = Offpsring()
%                 
%             Population(P(g_old>=g_new),:) = Offspring;
        end
        
        Archive = [Archive;[Population,FunctionValue]];
        
%         if mod(size(Archive,1),N)~=0
%             error('size of the archive must be a factor of N');
%         end
%          Solutions_gen(Gene+1).data = [Population,FunctionValue];
    end
end