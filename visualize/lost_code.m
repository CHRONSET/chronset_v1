% %%
% figure;
% subplot(131);
% a = gca;
% plot(sort(dum(:,2)-dum(:,1)),'b.');
% subplot(132);
% a = [a gca];
% plot(sort(dum(:,3)-dum(:,1)),'b.');
% subplot(133);
% a = [a gca];
% plot(sort(dum(:,3)-dum(:,2)),'b.');
% 
% lm1  = zeros(length(a),2);
% for it = 1:length(a)
%     hold(a(it),'on');
%     plot(a(it),[1 length(dum(:,1))],[-10 -10],'r-');
%     plot(a(it),[1 length(dum(:,1))],[10 10],'r-');
%     lm1(it,:) = [get(a(it),'YLim')];
% end;
% axis(a,'tight');
% set(a,'YLim',[min(lm1(:,1)) max(lm1(:,2))]);
% 
% [d1] = dum(:,1)-dum(:,2);
% [d2] = dum(:,1)-dum(:,3);
% [d3] = dum(:,2)-dum(:,3);
% 
% [sel_idx] = intersect(intersect(find(d1 >= -50 & d1 <= 50),find(d2 >= -50 & d2 <= 50)),find(d3 >= -50 & d3 <= 50));
% 
% on = on(sel_idx);
% dum = dum(sel_idx,:);
% data = data(sel_idx,:);
% %%
% figure;
% subplot(131);
% a = gca;
% plot(dum(:,5),dum(:,1),'b.');
% subplot(132);
% a = [ gca];
% plot(dum(:,5),dum(:,2),'b.');
% subplot(133);
% a = [ gca];
% plot(dum(:,5),dum(:,4),'b.');
% 
% for it = 1:length(a)
%     xlabel(a(it),'SayWhen [ms]');
%     ylabel(a(it),'Manual [ms]');
% end;
% axis(a,'tight');
% set(gcf,'Color','w');
% %% 
% figure;
% 
% subplot(1,3,1);
% 
% plot(diag(S)./sum(diag(S)),'bs-','MarkerFaceColor','b');
% axis tight;
% ylim([0 1]);
% xlabel('Factor');
% ylabel('% of variance explained');
% set(gca,'XTick',[1:5]);
% set(gca,'YTick',[0 .5 1]);
% xlim([0 size(V,2)+1]);
% 
% mY = median(pred,2);
% 
% X = dum(:,[1,2,4]);
% MSb = zeros(1,size(X,1));
% for kt = 1:size(X,1)
%     MSb(kt) = sum((mean(X,2)-mean(X(kt,:),2)).^2)/(size(X,1)-1);
% end;
% 
% MSw = sum(sum((X-repmat(mean(X,2),[1 size(X,2)])).^2,2),1)/(size(X,1)*(size(X,2)-1));
% 
% ICC = (MSb-MSw)/(MSb+(size(X,2)-1)*MSw);
% d1 = [X(:,1)-X(:,2) X(:,2)-X(:,3) X(:,1)-X(:,3)];
% 
% X = pred;
% MSb = zeros(1,size(X,1));
% for kt = 1:size(X,1)
%     MSb(kt) = sum((mean(X,2)-mean(X(kt,:),2)).^2)/(size(X,1)-1);
% end;
% 
% MSw = sum(sum((X-repmat(mean(X,2),[1 size(X,2)])).^2,2),1)/(size(X,1)*(size(X,2)-1));
% 
% ICC2 = (MSb-MSw)/(MSb+(size(X,2)-1)*MSw);
% d2 = [X(:,1)-X(:,2) X(:,2)-X(:,3) X(:,1)-X(:,3)];
% 
% subplot(1,3,2);
% 
% bar([1 2],[ICC ICC2],'FaceColor',[.75 .75 .75],'EdgeColor','k');    
% ylim([.996 1]);
% xlim([0 3]);
% ylabel('Inter class correlation (ICC)');
% set(gca,'XTick',1:2);
% set(gca,'XTickLabel',{'Raw' 'SVD'});
% 
%     
% subplot(1,3,3);
% hold on;
% plot([1 length(d1(:))],[-50 -50],'--','Color',[.75 .75 .75]);
% plot([1 length(d1(:))],[50 50],'--','Color',[.75 .75 .75]);
% 
% plot(sort(d1(:)),'bs');    
% plot(sort(d2(:)),'ro');
% axis tight;
% 
% set(gcf,'Color','w');
% %%
% [U,S,V] = svd(dum(:,[1 2 3]));
% 
% recon = U(:,1)*S(1,1)*V(:,1)';
% tMe = mean(recon,2);
% 
% %%
% [b1,~,r1,~,~] = regress(tMe,[ones(size(dum(:,1))) dum(:,1)]);
% [b2,~,r2,~,~] = regress(tMe,[ones(size(dum(:,1))) dum(:,2)]);
% [b3,~,r3,~,~] = regress(tMe,[ones(size(dum(:,1))) dum(:,3)]);
% 
% yp1 = b1(1)+b1(2).*dum(:,1);
% yp2 = b2(1)+b2(2).*dum(:,2);
% yp3 = b3(1)+b3(2).*dum(:,3);
% 
% i1 = find(r1 >= -500 & r1 <= 500);
% i2 = find(r2 >= -500 & r2 <= 500);
% i3 = find(r3 >= -500 & r3 <= 500);
% 
% sel_idx = intersect(intersect(intersect(i1,i2),intersect(i1,i3)),intersect(i2,i3));
% 
% figure;
% subplot(131);
% a = gca;
% subplot(132);
% a = [a gca];
% subplot(133);
% a = [a gca];
% 
% for it = 1:length(a)
%     hold(a(it),'on');
%     plot(a(it),[1 max(max([yp1 yp2 yp3]))],[0 0],'r-');
% end;
% 
% plot(a(1),yp1,r1,'b.');
% plot(a(2),yp2,r2,'b.');
% plot(a(3),yp3,r3,'b.');
% 
% plot(a(1),yp1(sel_idx),r1(sel_idx),'ro');
% plot(a(2),yp2(sel_idx),r2(sel_idx),'ro');
% plot(a(3),yp3(sel_idx),r3(sel_idx),'ro');
% 
% axis(a,'tight');
% lm = get(a,'YLim');
% 
% set(a,'YLim',[min(sort([lm{:}])) max(sort([lm{:}]))]);
% 
% dum2 = [dum tMe];
% 
% X = dum2(setdiff(1:length(r1),sel_idx),:);
% MSb = zeros(1,size(X,1));
% for kt = 1:size(X,1)
%     MSb(kt) = sum((mean(X,2)-mean(X(kt,:),2)).^2)/(size(X,1)-1);
% end;
% 
% MSw = sum(sum((X-repmat(mean(X,2),[1 size(X,2)])).^2,2),1)/(size(X,1)*(size(X,2)-1));
% 
% ICC = (MSb-MSw)/(MSb+(size(X,2)-1)*MSw);
% 
% X = dum2(sel_idx,:);
% MSb = zeros(1,size(X,1));
% for kt = 1:size(X,1)
%     MSb(kt) = sum((mean(X,2)-mean(X(kt,:),2)).^2)/(size(X,1)-1);
% end;
% 
% MSw = sum(sum((X-repmat(mean(X,2),[1 size(X,2)])).^2,2),1)/(size(X,1)*(size(X,2)-1));
% 
% ICC2 = (MSb-MSw)/(MSb+(size(X,2)-1)*MSw);