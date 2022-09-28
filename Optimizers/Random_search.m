function [Population,FunctionValue,Archive]=Random_search(Problem,MaxFun,num_dimensions,Population,N)
    
% Solutions_gen = struct;
addpath(genpath('Public'));
M = Problem.numberOfObjectives;

MaxValue = 1*ones(1,num_dimensions);
MinValue = -1*ones(1,num_dimensions);
    
FunctionValue = zeros(size(Population,1),M);
for i = 1:size(Population,1)
    FunctionValue(i,:) = Problem.evaluate(Population(i,:)); 
end
Archive = [Population,FunctionValue]; % storing all evaluated solutions
    
% Solutions_gen(1).data = [Population,FunctionValue]; 
%%  
nfun = 0;
    %% Optimization
 while nfun <= MaxFun    

    Offspring = repmat(MinValue,N,1) + (repmat(MaxValue,N,1) - repmat(MinValue,N,1)).*rand(N,num_dimensions);
    Fitness = zeros(size(Offspring,1),M);

    for i = 1:size(Offspring,1)
        Fitness(i,:) = Problem.evaluate(Offspring(i,:));
    end
%     Archive = [Archive;[Offspring,Fitness]];

    nfun = nfun + size(Offspring,1);

Population = Offspring;
FunctionValue = Fitness;
% Solutions_gen(Gene+1).data = [Population,FunctionValue];
Archive = [Archive;[Population,FunctionValue]];

% if mod(size(Archive,1),N)~=0
%             error('size of the archive must be a factor of N');
% end

end

end
