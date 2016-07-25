function [params] = compute_slope_parameters4onsets(SA,phID,n)

%%
rsq = zeros(length(SA),1);
betas = zeros(length(SA),2);
int_b = zeros(length(SA),4);

yp = cell(length(SA),1);
for it = 1:length(SA)

    %SA{it} = SA{it}(find(sign(SA{it}(:,1)-250)==1),:);
    X = SA{it}(:,1);
    Y = SA{it}(:,2);
    % todo: adjust CI for family-wise comps 1-p^n
    fwa = 1-(1-(.95^((1/(length(phID)*2)))));
    [betas(it,:),bint,~,~,stats]= regress(Y,[ones(size(X)) X],1-fwa);
    int_b(it,:) = bint(:)';
    
    yp{it} = betas(it,1)+betas(it,2).*X;
    rsq(it) = stats(1);
end;
params.SA = SA;
params.yp = yp;
params.rsq = rsq;
params.betas = betas;
params.int_b = int_b;
params.phID = phID;
params.n = n;