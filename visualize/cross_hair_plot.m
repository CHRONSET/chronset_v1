function cross_hair_plot(params)

%%
chck = zeros(size(params.int_b,1),2);

for it = 1:size(params.int_b,1);
   chck(it,1) = sum(sign(params.int_b(it,[2 4])-mean(params.betas(:,2))));
   chck(it,2) = sum(sign(params.int_b(it,[1 3])-mean(params.betas(:,1))));
end;
hold on;
%plot(params.betas(:,1),params.yp,'r','LineWidth',3);
%plot([min(params.int_b(:,1)) max(params.int_b(:,3))],[1 1],'b','LineWidth',3);
%plot([0 0],[min(params.int_b(:,2)) max(params.int_b(:,4))],'b','LineWidth',3);
%%
%plot([min(params.int_b(:,1)) max(params.int_b(:,3))],repmat(mean(params.betas(:,2)),[1 2]),'c','LineWidth',3);
%plot(repmat(mean(params.betas(:,1)),[1 2]),[min(params.int_b(:,2)) max(params.int_b(:,4))],'c','LineWidth',3);

[b,~,~,~,~] = regress(params.betas(:,2),[ones(size(params.betas(:,2))) params.betas(:,1)]);
yp = b(1)+b(2).*params.betas(:,1);
for it = 1:size(params.betas,1)
    x = params.int_b(it,:);
    x = reshape(x,[2 2]);
    if ismember(it,find(min(chck,[],2)==0))
        plot(repmat(mean(x(1,:)),[1 2]),[x(2,1) x(2,2)],'Color',[.75 .75 .75]);
        plot([x(1,1) x(1,2)],repmat(mean(x(2,:)),[1 2]),'Color',[.75 .75 .75]);
        plot(params.betas(it,1),params.betas(it,2),'k.');
        %text(params.betas(it,1)+2,params.betas(it,2)+.02,params.phID{it},'Color',[.75 .75 .75]);
    end;
end;

for it = 1:size(params.betas,1)
    x = params.int_b(it,:);
    x = reshape(x,[2 2]);
    if ismember(it,find(min(chck,[],2)==0))
    else
        plot(repmat(mean(x(1,:)),[1 2]),[x(2,1) x(2,2)],'k-','LineWidth',3);
        plot([x(1,1) x(1,2)],repmat(mean(x(2,:)),[1 2]),'k-','LineWidth',3);
        plot(params.betas(it,1),params.betas(it,2),'.','Color',[.75 .75 .75]);
    end;
end;


for it = 1:size(params.betas,1)    
    if ismember(it,find(min(chck,[],2)==0))
    else
        text(params.betas(it,1)+2,params.betas(it,2)+.02,params.phID{it},'FontWeight','bold','Color','r');

    end;
end;

plot([mean(params.betas(:,1))-std(params.betas(:,1)) mean(params.betas(:,1))+std(params.betas(:,1))],repmat(mean(params.betas(:,2)),[1 2]),'c','LineWidth',3);
plot(repmat(mean(params.betas(:,1)),[1 2]),[mean(params.betas(:,2))-std(params.betas(:,2)) mean(params.betas(:,2))+std(params.betas(:,2))],'c','LineWidth',3);


axis tight;
%xlim([-70 max(max(params.int_b(:,[1 3])))]);
%ylim([.85 1.15]);
xlabel('Intercept');
ylabel('Slope');