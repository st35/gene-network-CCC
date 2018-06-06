% This code calculates the normalized CCC for samples from 2 phenotypic groups.

clc;
clear all;

data = importdata('expression_data.txt', ',');
% Import gene expression data for different samples. The file expression_data.txt must be CSV file with gene names along rows and sample names along columns.
% Only genes for which you want to calculate the CCC must be included in the expression_data.txt file.
expression_data = data.data;
expression_data = transpose(expression_data);
expression_data = log2(expression_data); % Do this only if the expression data is not log-normalized.
      
num_genes = size(expression_data, 2); % Get number of genes in the dataset.
num_pairs = nchoosek(num_genes, 2); % Get number of gene pairs.

grp_1_indices = importdata('group_1_indices.txt'); % Column indices of samples belonging to phenotypic group 1 in the file expression_data.txt.
grp_2_indices = importdata('group_2_indices.txt'); % Column indices of samples belonging to phenotypic group 2 in the file expression_data.txt.

num_bootstrap_samples = 100; % Number of bootstrap samples you want to draw.

% CCC calculation for phenotypic group 1.
grp_1_expression = expression_data(grp_1_indices, :);
grp_1_size = size(grp_1_expression, 1);

% Bootstrap method to obtain the error in the estimate of normalized CCC.
CCC_grp1 = zeros(1, num_bootstrap_samples);

for T = 1:num_bootstrap_samples
    bootstrap_samp = randsample(1:grp_1_size, grp_1_size, true); % Draw samples randomly with replacememt/
    bootstrap_samp_data = grp_1_expression(bootstrap_samp, :);
    CCC = get_CCC(bootstrap_samp_data); % Get CCC of the bootstrap sample drawn.
    CCC_rand = get_CCC_random_net(bootstrap_samp_data, 10); % Get mean CCC of 10 random nets with the same set of nodes but with randomly re-distributed edge weights.
    CCC_grp1(1, T) = (CCC - CCC_rand) / (1 - CCC_rand); % Get normalized CCC of the bootstrap sample.
end

full_samp = 1:grp_1_size; % To calculate the actual CCC for phenotypic group 1.
full_samp_data = grp_1_expression(full_samp, :);
CCC = get_CCC(full_samp_data);
CCC_rand = get_CCC_random_net(full_samp_data, 10);
final_CCC_grp1 = (CCC - CCC_rand) / (1 - CCC_rand);

% CCC calculation for phenotypic group 2.
grp_2_expression = expression_data(grp_2_indices, :);
grp_2_size = size(grp_2_expression, 1);

% Bootstrap method to obtain the error in the estimate of normalized CCC.
CCC_grp2 = zeros(1, num_bootstrap_samples);

for T = 1:num_bootstrap_samples
    bootstrap_samp = randsample(1:grp_2_size, grp_2_size, true); % Draw samples randomly with replacememt/
    bootstrap_samp_data = grp_2_expression(bootstrap_samp, :);
    CCC = get_CCC(bootstrap_samp_data); % Get CCC of the bootstrap sample drawn.
    CCC_rand = get_CCC_random_net(bootstrap_samp_data, 10); % Get mean CCC of 10 random nets with the same set of nodes but with randomly re-distributed edge weights.
    CCC_grp2(1, T) = (CCC - CCC_rand) / (1 - CCC_rand); % Get normalized CCC of the bootstrap sample.
end

full_samp = 1:grp_2_size; % To calculate the actual CCC for phenotypic group 2.
full_samp_data = grp_2_expression(full_samp, :);
CCC = get_CCC(full_samp_data);
CCC_rand = get_CCC_random_net(full_samp_data, 10);
final_CCC_grp2 = (CCC - CCC_rand) / (1 - CCC_rand);

% Plot the results.
CCC_estimate = [final_CCC_grp1 final_CCC_grp2];
error_estimate = [std(CCC_grp1) std(CCC_grp2)];

figure;
hold on;
bar(1:2, CCC_estimate, 'w');
errorbar(1:2, CCC_estimate, error_estimate, '.');
xticks(1:2);
set(gca, 'xticklabel', {'Group 1', 'Group 2'});
ylabel('Normalized CCC');
set(gca, 'FontSize', 16);

% Estimate the p-value. Here, the null hypothesis is that the CCC of group 1 is not higher than the CCC of group 1.

d1 = CCC_grp1;
d2 = CCC_grp2;

u = final_CCC_grp1 - final_CCC_grp2;

z = (final_CCC_grp1 + final_CCC_grp2) / 2;

a1 = d1 - final_CCC_grp1 + z;
a2 = d2 - final_CCC_grp2 + z;

ub = a1 - a2;

p_value = sum(ub > u) / length(d1);