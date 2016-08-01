%%
path2files = '~/froux/chronset/data/BCN/feature_data/';
files = dir([path2files,'*.mat']);

x1 =zeros(length(files),6);
r1 =zeros(length(files),6,6);

parfor jt = 1:length(files)
    
    dat = load([path2files,files(jt).name]);
    
    M = ([dat.savedata.features{:}]);

    [U,S,V] =svd((M'*M));
    
    x1(jt,:) = cumsum(diag(S))./sum(diag(S));
    r1(jt,:,:) = corr(M);

end;
%%
path2files = '~/froux/chronset/data/SayWhen/feature_data/';
files = dir([path2files,'*.mat']);

x2 =zeros(length(files),6);
r2 =zeros(length(files),6,6);

parfor jt = 1:length(files)
    
    dat = load([path2files,files(jt).name]);
    
    M = ([dat.savedata.features{:}]);

    [U,S,V] =svd((M'*M));
    
    x2(jt,:) = cumsum(diag(S))./sum(diag(S));
    r2(jt,:,:) = corr(M);
    
end;
%%
figure;
subplot(221);
imagesc(squeeze(mean(r1,1)));
set(gca,'XTick',[1:6]);
set(gca,'YTick',[1:6]);
caxis([-1 1]);

subplot(222);
errorbar(1:6,mean(x1,1),std(x1,0,1),'bs-','MarkerFaceColor','b');ylim([min(mean(x1,1)-std(x1,0,1))-.01 1.02]);
set(gca,'XTick',[1:6]);
xlim([0 7]);

subplot(223);
imagesc(squeeze(mean(r2,1)));
set(gca,'XTick',[1:6]);
set(gca,'YTick',[1:6]);
caxis([-1 1]);

subplot(224);
errorbar(1:6,mean(x2,1),std(x2,0,1),'bs-','MarkerFaceColor','b');ylim([min(mean(x2,1)-std(x2,0,1))-.01 1.02]);
set(gca,'XTick',[1:6]);
xlim([0 7]);



