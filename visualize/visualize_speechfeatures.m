function visualize_speechfeatures(feat_data,features,finf)
%%
%load('/bcbl/home/home_a-f/froux/BilingualNamingExp/results/Colormap_spectrogram.mat');
%%
feat_data.t = feat_data.t*1000;
feat_data.tAx = feat_data.tAx.*1000;
figure;
subplot(3,1,1:2);
a = gca;
%imagesc(round(feat_data.t),feat_data.f(2:end),squeeze(mean(feat_data.ds,1))');
Y = squeeze(feat_data.ds(2,:,:))';
% Y = Y-min(min(Y));
% Y = Y./(max(max(Y))-min(min(Y)));
imagesc(round(feat_data.t),feat_data.f(2:end),Y);
axis xy;
caxis([-.1e-6 .1e-6]);
colormap gray;
%colormap(c);
hold on;
set(gca,'XTick',(round([0:0.250:max(feat_data.t)].*1000)));

a(2) = axes('position',[.13 .51 .776 .2]);
hold on;
h = zeros(1,length(finf.xidx_p));
for it = 1:length(finf.xidx_p)
    x = feat_data.t;
    y = features{finf.inf_f(1)}.*(finf.x_tresh{finf.inf_f(1)})';
    y(y==0) = NaN;
    plot(x,y,'r-','LineWidth',4);
end;
% k = 0;
% for it = 1:length(finf.xidx_p2)
%     k = k+1;
%     h(k) = area(feat_data.t(finf.xidx_p2(it):finf.xidx_n2(it)),ones(1,length(finf.xidx_p2(it):finf.xidx_n2(it))),-2);
% end;
% set(h,'FaceColor','g');
% set(h,'EdgeColor','w');
ylim([min(y) max(y)]);
%
plot([feat_data.t(1) feat_data.t(end)],[finf.tresh{finf.inf_f(1)} finf.tresh{finf.inf_f(1)}],'m--','LineWidth',6)
xlim([feat_data.t(1) feat_data.t(end)]);
ylim([min(y) max(y)+.1]);
axis off;

a(3) = axes('position',[.13 .72 .776 .2]);
hold on;
for it = 1:length(finf.xidx_n)
    x = feat_data.t;
    y = features{finf.inf_f(2)}.*(finf.x_tresh{finf.inf_f(2)})';
    y(y==0) = NaN;
    plot(x,y,'y-','LineWidth',4);
end;
plot([feat_data.t(1) feat_data.t(end)],[finf.tresh{finf.inf_f(2)} finf.tresh{finf.inf_f(2)}],'g--','LineWidth',6);
axis tight;
xlim([feat_data.t(1) feat_data.t(end)]);
axis off;

a(4) = axes('position',[.13 .429 .776 .2]);
hold on;
for it = 1:length(finf.xidx_n)
    x = feat_data.t;
    y = features{finf.inf_f(3)}.*(finf.x_tresh{finf.inf_f(3)})';
    y(y==0) = NaN;
    plot(x,y,'c-','LineWidth',4);
end;
plot([feat_data.t(1) feat_data.t(end)],[finf.tresh{finf.inf_f(3)} finf.tresh{finf.inf_f(3)}],'b--','LineWidth',6);
axis tight;
xlim([feat_data.t(1) feat_data.t(end)]);
axis off;
%%
subplot(3,1,3);
hold on;
plot(feat_data.tAx,feat_data.orig,'k');
plot([feat_data.t(finf.on_idx) feat_data.t(finf.on_idx)],[min(feat_data.orig) max(feat_data.orig)],'r','LineWidth',4);
plot([feat_data.t(finf.off_idx) feat_data.t(finf.off_idx)],[min(feat_data.orig) max(feat_data.orig)],'r','LineWidth',4);
%plot([feat_data.t(finf.on_idx) feat_data.t(finf.off_idx)],[min(feat_data.orig) max(feat_data.orig)],'r','LineWidth',4);
%plot([feat_data.t(finf.on_idx) feat_data.t(finf.off_idx)],[max(feat_data.orig) min(feat_data.orig)],'r','LineWidth',4);
set(gca,'XTick',(round([0:0.250:max(feat_data.t)].*1000)));
% h = zeros(1,length(finf.xidx_p2));
% k = 0;
% for it = 1:length(finf.xidx_p2)
%     k = k+1;
%     h(k) = area(feat_data.t(finf.xidx_p2(it):finf.xidx_n2(it)),(min(feat_data.orig)+.5).*ones(1,length(finf.xidx_p2(it):finf.xidx_n2(it))),min(feat_data.orig));
% end;
% set(h,'FaceColor','g');
% set(h,'EdgeColor','k');
ylim([min(feat_data.orig) max(feat_data.orig)]);
set(gca,'Fontsize',20);
axis tight;
ylabel(gca,'Speech amplitude [\sigma]','Fontsize',20);
xlabel(gca,'Time [ms]','Fontsize',20);

%set(gca,'XTick',(round([min(feat_data.t):0.25:max(feat_data.t)].*100))/100);
% %%
% subplot(3,1,3);
% [x,y] = ginput(2);
% idx = find(feat_data.tAx >=x(1) & feat_data.tAx <=x(2));
%
% p1 = audioplayer(feat_data.orig(idx),feat_data.FS);
% play(p1);
set(gcf,'Color','w');
set(a,'Fontsize',20);
yl = get(a,'YTick');
yl = {yl{1}./10e2};
set(a(1),'YTickLabel',yl);
ylabel(a(1),'Frequency [KHz]','Fontsize',20);
xlabel(a(1),'Time [ms]','Fontsize',20);