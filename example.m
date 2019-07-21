clear
load example.mat % Load a Author-Conference dataset
Net = C_Net; % Load network
P = Net * diag(sum(Net))^-1; % Get transition matrix
%%%
% [1, 20] are two seed nodes 1 and 20
% [1, 2] means the seeds nodes 1 and 20 are of different colors
% 0.9 is the forwarding probability, which is always 0.9 in my program
% @(t)1 is an anonymous decaying function, which always return 1
% 10 is the number of iterations, which is always 10
% 10000 and -10000 are attraction and repulsion strength
R = crw_sparse( P, [1, 20], [1, 2], 0.9, @(t)1, 10, 100, -100);
R = reshape(R, 20, 2);
S1 = sweep_mex(P, R(:, 1), 5); % Ground Truth [1; 2; 7; 8; 15], 80% hit
S2 = sweep_mex(P, R(:, 2), 5); % Ground Truth [6; 9; 11; 19; 20], 100% hit