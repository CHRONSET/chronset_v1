%%
function [n,x] = compute_pct_trial(D,bs)
[n,x] = hist(D,bs);
n = n./sum(n);