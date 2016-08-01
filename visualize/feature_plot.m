visualize_speechfeatures(feat_data,feat_data.features,feat_data.finf);
%%
set(gcf,'PaperPositionMode','auto');
savepath = '';
print(gcf,'-dtiff','-r300','-cmyk',[savepath,'spectrogram_1.tiff']);
%%
figure;
plot_speech_features(feat_data);
%%
set(gcf,'PaperPositionMode','auto');
savepath = '';
print(gcf,'-dtiff','-r1200','-zbuffer',[savepath,'features.tiff']);