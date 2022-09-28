function Sets = generate_sets_Indicator(Data,no_sol,no_var,ref_point)

F = Data(:,no_var+1:end);

obj = F;
N = size(obj,1);
PopObj = (obj-repmat(min(obj),N,1))./(repmat(max(obj)-min(obj),N,1));
I = zeros(N);
for i = 1 : N
    for j = 1 : N
        I(i,j) = max(PopObj(i,:)-PopObj(j,:));
    end
end
C = max(abs(I));
Fitness = sum(-exp(-I./repmat(C,N,1)/kap)) + 1;

[~,id] = sort(Fitness);
id_store = id;

j = 1;
while ~isempty(Data)
    
    if size(Data,1)<no_sol
        rem = no_sol - size(Data,1);
        index1 = id(1:no_sol);
        id_2 = index1(end:end-length(index1));
        index2 = id_store(id_2:id_2+rem);
        Sets{j}.X = Data(index2,1:no_var);
        Sets{j}.F = Data(index2,no_var+1:end);
        Sets{j}.HV = Fitness(index2,:);
    else
        index = id(1:no_sol);
        Sets{j}.X = Data(index,1:no_var);
        Sets{j}.F = Data(index,no_var+1:end);
        Sets{j}.HV = Fitness(index,:);
    
        Data(index,:) = [];
    
        id(index,:) = [];
    end
    j = j+1;
end

