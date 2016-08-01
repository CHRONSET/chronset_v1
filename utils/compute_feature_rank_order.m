function [n,x] = compute_feature_rank_order(inf_f)

%%
n = zeros(6,6);
x = zeros(6,6);
for it = 1:6%loop over rankings
    
    % count feature code per ranking
    [n(it,:),] = hist(inf_f(:,it),[1:6]);
    
end;

%n has dims ranking x features and contains count values

% normalize count values to % trials
for it = 1:size(n,1)
    
    n(it,:) = n(it,:)./sum(n(it,:));
end;

x = 1:size(inf_f,2);

