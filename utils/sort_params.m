function [params] = sort_params(params,idx)

params.SA = params.SA(idx);
params.yp = params.yp(idx);
params.rsq = params.rsq(idx);
params.betas = params.betas(idx,:);
params.int_b = params.int_b(idx,:);
params.phID = params.phID(idx,:);
params.n = params.n(idx,:);