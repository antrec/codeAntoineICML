function [simMatrix] = gen_sim_from_ecoli()
% generate similarity matrix with a few outliers from E. coli DNA data

X = csvread('ecoli_coo_mat.csv');
n = max(X(:,1));
simMatrix = sparse(X(:,1), X(:,2), X(:,3), n, n);