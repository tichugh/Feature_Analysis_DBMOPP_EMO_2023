%% script for running different MOEAs on DBMOPP instances
%% The instances are in 'lille_experiment_instances.mat'

% Please contact t.chugh@exeter.ac.uk if you find any bug or problem in
% running the code


clear;
clc; close all;
% Optimizer = 'NSGAII';
% Optimizer = 'IBEA';
Optimizer = 'Random';
lhs_dir = 'search-space';
addpath(genpath('Optimizers'));
Results_folder = ([Optimizer '_Results']);


load lille_experiment_instances.mat;


%%
MaxFun = 50000;
MaxRun_optimizers =30; % number of folds per instance
pop_sizes = [100, 105, 120, 126, 132, 112, 156, 90, 275]; % population size - depends on number of objectives

%% 
for i  = 1:length(instances)
    for run = 1:MaxRun_optimizers
        Problem = instances{i,run};
        if ~isempty(Problem)
            n_var = Problem.numberOfDesignVariables;
            n_obj = Problem.numberOfObjectives;
            pop_size = pop_sizes(n_obj-1);
            size_lhs = 200*n_var;
            data_lhs = readtable([lhs_dir '\lhs_d' num2str(n_var) '_n' num2str(size_lhs) '.csv']);
            Initial_Population = table2array(data_lhs);       
            if strcmp(Optimizer,'NSGAII')
                [~,~,Archive] = NSGA_II(Problem,MaxFun,n_var,Initial_Population,pop_size);
                csvwrite([Results_folder '\Solutions_NSGAII_id_' num2str(i) '_run_' num2str(run) '.csv'],Archive);
            elseif strcmp(Optimizer,'MOEAD')
                [~,~,Archive] = MOEAD(Problem,MaxFun,n_var,Initial_Population,pop_size);
                csvwrite([Results_folder '\Solutions_MOEAD_id_' num2str(i) '_run_' num2str(run) '.csv'],Archive);
            elseif strcmp(Optimizer,'IBEA')
                [~,~,Archive] = IBEA(Problem,MaxFun,n_var,Initial_Population,pop_size);
                 csvwrite([Results_folder '\Solutions_IBEA_id_' num2str(i) '_run_' num2str(run) '.csv'],Archive);
            elseif strcmp(Optimizer,'Random')
                [~,~,Archive] = Random_search(Problem,MaxFun,n_var,Initial_Population,pop_size);
                csvwrite([Results_folder '\Solutions_Random_id_' num2str(i) '_run_' num2str(run) '.csv'],Archive);
            end
        end
    end
end