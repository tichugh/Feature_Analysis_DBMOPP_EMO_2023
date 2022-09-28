function [Population,FunctionValue,Archive]=IBEA(Problem,MaxFun,num_dimensions,Population,N)
    
% Solutions_gen = struct;
addpath(genpath('IBEA_files'));
addpath(genpath('Public'));
%     M = Problem.num_objectives;
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
    FrontValue = P_sort(FunctionValue);
    CrowdDistance = F_distance(FunctionValue,FrontValue); % crowding distance calculation for NSGA-II
 %%   
    nFun = 0;
    while nFun <= MaxFun     
%         Gene
%         MatingPool = F_mating(Population,FrontValue,CrowdDistance); % mating pool
        MatingPool = F_mating(Population,N); % mating pool
        Coding = 'Real';
        Offspring = P_generator(MatingPool,Boundary,Coding,N);
        Fitness = zeros(size(Offspring,1),M);
    
        for i = 1:size(Offspring,1)
            Fitness(i,:) = Problem.evaluate(Offspring(i,:));
        end
        
        nFun = nFun + size(Offspring,1);
%         Archive = [Archive;[Offspring,Fitness]];
        
        Population = [Population;Offspring];
        FunctionValue = [FunctionValue;Fitness];

        [Population,FunctionValue]=IBEA_Selection(Population,N,0.5,FunctionValue);
%         Solutions_gen(Gene+1).data = [Population,FunctionValue];
        Archive = [Archive;[Population,FunctionValue]];
        
%         if mod(size(Archive,1),N)~=0
%             error('size of the archive must be a factor of N');
%         end
    end

end
