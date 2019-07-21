function [ Rt ] = crw_sparse( P, Seeds, Colors, alpha, psi, n_iters, L1, L2)
% crw : Power iteration version of colored rw
%%% Input
% P: transition matrix, col-stochastic
% Seeds: cell array of seed nodes with constrains
% alpha: forward probability
% n_iters: max iteration numbers
% L1: factor of same color, in (1, infinity)
% L2: factor of diff color, in (-1, -0)
%%% Output
% R: Color score matrix
K = length(unique(Colors));
N = size(P, 1);
% Prepare s
S = sparse(zeros(N*K, 1));
for k = 1:K
    seed = Seeds(Colors==k);
    S((k-1)*N+seed) = 1/numel(seed);
end
Rt = S;
P0 = kron(eye(K), P);
PP = P0; % PP === \mathbb{P}
% Lambda = sparse(diag(repmat(L1, 1, K)));
Lambda = ones(K, K) * L2;
Lambda = Lambda - eye(K) * L2 + eye(K) * L1;
Lambda = kron(Lambda, speye(N));
for t = 1:n_iters
    Pt = Lambda * Rt;
    Pt = reinforce(P0, full(Pt));
    x = psi(t);
    PP = x * Pt + (1-x) * PP;
    Rt = alpha * PP * Rt + (1-alpha) * S;
end
end