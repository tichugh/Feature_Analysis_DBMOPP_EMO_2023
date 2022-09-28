function [Population,FunctionValue] = IBEA_Selection(Population,N,kap,FunctionValue)
Next = 1 : length(Population);
[Fitness,I,C] = Fitness_IBEA(FunctionValue,kap);
    while length(Next) > N
        [~,x]   = min(Fitness(Next));
        Fitness = Fitness + exp(-I(Next(x),:)/C(Next(x))/kap);
        Next(x) = [];
    end
    Population = Population(Next,:);
    FunctionValue = FunctionValue(Next,:);
end

function [Fitness,I,C] = Fitness_IBEA(obj,kap)
N = size(obj,1);
    PopObj = (obj-repmat(min(obj),N,1))./(repmat(max(obj)-min(obj),N,1));
    I      = zeros(N);
    for i = 1 : N
        for j = 1 : N
            I(i,j) = max(PopObj(i,:)-PopObj(j,:));
        end
    end
    C = max(abs(I));
    Fitness = sum(-exp(-I./repmat(C,N,1)/kap)) + 1;
    
end