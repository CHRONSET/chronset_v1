%%
load('/bcbl/home/home_a-f/froux/chronset/thresholds/mc_confidence_estimates_01-Oct-2015.mat');

figure;
a = zeros(size(T1,1),1);
h = zeros(size(T1,2),1);
for it = 1:size(T1,1)
    
    subplot(1,2,it);
    a(it) = gca;
    hold on;
    plot([5 50],[0 0],'k--');
    
    for jt = 1:size(T1,2)
        
        [x1] = squeeze(T1(it,jt,:));
        [x2] = squeeze(T2(it,jt,:));
        
        
        if jt == 1
            c = {[.9 0 0],'-o'};
        elseif jt == 2
            c = {[0 .75 .9],'-s'};
        elseif jt == 3
            c = {[.75 .75 .75],'-^'};
        else
            c = {[.75 0 .75],'-^'};
        end;
        
        h(jt) = plot(5:50,(x1-x2)./x1,[c{2}],'Color',[c{1}],'MarkerFaceColor',[c{1}],'MarkerSize',4);
        
    end;
    
end;
axis(a,'tight');
set(a,'YLim',[-.005 .6]);
for it = 1:length(a)
    xlabel(a(it),'Sample size');
    ylabel(a(it),'Reduction in statistical power [%]');
end;
legend(a(1),h,'SD:50','SD:75','SD:100','SD:125');
title(a(1),'Effect size: 15 ms');
title(a(2),'Effect size: 30 ms');
set(gcf,'Color','w');

set(a,'LineWidth',3);
