function [ ids ] = spectral_clustering(W, k)

D = diag(1./(eps+sqrt(sum(W, 2))));
W = D * W * D;
[U, s, V] = svd(W);
V = U(:, 1 : k);
V = normr(V);
ids=litekmeans(V, k,  'Replicates', 500);
end
