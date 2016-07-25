function lattice_plot_1(params)

[~,s_idx] = sort(params.rsq);

params.SA = params.SA(s_idx);
params.yp = params.yp(s_idx);
params.rsq = params.rsq(s_idx);
params.betas = params.betas(s_idx,:);
params.int_b = params.int_b(s_idx,:);
params.phID = params.phID(s_idx);
params.n = params.n(s_idx);
%%
for it = 1:length(params.SA)
    params.SA{it} = params.SA{it}./1000;
    params.yp{it} = params.yp{it}./1000;
    
    x = params.SA{it}(:,1);
    x = (x-min(x))./(max(x)-min(x));
    params.SA{it}(:,1) = x;
    
    x = params.SA{it}(:,2);
    x = (x-min(x))./(max(x)-min(x));
    params.SA{it}(:,2) = x;
    
    x = params.yp{it};
    x = (x-min(x))./(max(x)-min(x));
    params.yp{it} = x;    
end;
%%
lm = zeros(length(params.SA),2);
for it = 1:length(params.SA)    
    lm(it,:) = [min(min(params.SA{it})) max(max(params.SA{it}))];
end;
lm = [min(min(lm)) max(max(lm))];
lm = round(lm);
%%
figure;
subplot(round(length(params.SA)/5),5,1);
xlim([lm(1) lm(2)+.1]);
ylim([lm(1) lm(2)+.1]);
set(gca,'XTick',round([lm(1) lm(2)].*1)/1);
set(gca,'YTick',round([lm(1) lm(2)].*1)/1);

k = 1;
for it = 1:length(params.SA)
    
    k = k+1;
    subplot(round(length(params.SA)/5),5,k);
    hold on;
    plot(params.SA{it}(:,1),params.yp{it},'r');
    plot(params.SA{it}(:,1),params.SA{it}(:,2),'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
    axis tight;
    xlim([0 1+.1]);
    ylim([0 1+.1]);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    title([params.phID(it),'R^2:',round(params.rsq(it)*100)/100]);
end;
set(gcf,'Color','w');

return;