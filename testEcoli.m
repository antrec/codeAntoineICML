clear;
close all;
addpath(genpath(pwd));
rng(0);
A = gen_sim_from_ecoli();
n = size(A,1);

thisExpName = 'testAllMethodsEcoli.m';
% Chose parameter dh according to number of non-zero elements of A
nAll = floor(1/2 * nnz(A));
nInDiags = ((1:n)+1).*(2*n-(1:n))/2;
dh = find((nInDiags - nAll)>0, 1, 'first');

% Run all methods
[perms, huberscores, twosumscores, elTimes] = testAllMethods(A, dh);

% Save results
save(thisExpName, 'perms', 'huberscores', 'twosumscores', 'elTimes');