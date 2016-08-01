function plot_speech_features(feat_data)
%%
feat_data.t = feat_data.t.*1000;
feat_data.t = round(feat_data.t);
figure;
a = zeros(length(feat_data.features),1);
for it = 1:length(feat_data.features)
    
    subplot(3,2,it);
    a(it) = gca;
    hold on;
    plot(feat_data.t,feat_data.features{it},'k-','LineWidth',3);
    %plot([min(feat_data.t) max(feat_data.t)],[feat_data.finf.tresh{it} feat_data.finf.tresh{it}],'r--','LineWidth',3);
    axis tight;
    
end;

set(a,'Fontsize',20);
set(a,'YTick',[0 1]);
set(gcf,'Color','w');
set(a,'Box','on');

for it = 1:length(a)
    xlabel(a(it),'Time [ms]','Fontsize',20);
    ylabel(a(it),'[a.u]','Fontsize',20);
end;

title(a(1),'Amplitude','Fontsize',20);
title(a(2),'Wiener entropy','Fontsize',20);
title(a(3),'Spectral change','Fontsize',20);
title(a(4),'Amplitude modulation','Fontsize',20);
title(a(5),'Frequency modulation','Fontsize',20);
title(a(6),'Harmonic pitch','Fontsize',20);
