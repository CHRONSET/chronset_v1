%%
addpath('/bcbl/home/home_a-f/froux/mFiles/shadedErrorBar');
addpath('/bcbl/home/home_a-f/froux/BilingualNamingExp/mfiles/');
addpath('/bcbl/home/home_a-f/froux/BilingualNamingExp/mfiles/Grifith/');       
%% init workspace
clear;
clc;
%close all;
if matlabpool('size') ==0
    matlabpool;
end;
%% read in the manual scores
path2files = '/bcbl/home/public/YAZ/Expe_NP/NP/Mono_copy/semi_automatic/';
files = dir([path2files,'*.azk']);

RTdata = cell(length(files),1);
parfor it = 1:length(files)
    
    fid = fopen([path2files,files(it).name],'r+');
    [data] = textscan(fid,'%s');
    data = data{:};
    data = data(find(strcmp(data,'RT'))+1:end);
    
    data = str2double(data);
    
    [RTdata{it}] = reshape(data,[2 length(data)/2])';
    
end;
clear data;
%%
n = zeros(length(RTdata),1);
for it = 1:length(RTdata)
    n(it) = size(RTdata{it},1);
end;

path2files2 = '/bcbl/home/public/YAZ/Expe_NP/NP/Mono_copy/automatic_scores/';
dat = cell(sum(n),1);
aRT = cell(length(files),1);
for it = 1:length(files)
    
    dat = cell(n(it),1);
    parfor jt = 1:size(RTdata{it},1)
        files = dir([path2files,'*.azk']);
        f = dir([path2files2,files(it).name(1:regexp(files(it).name,'.azk')-1),'_*_',num2str(RTdata{it}(jt,1)),'.mat']);
        dum = load([path2files2,f.name]);
        dat{jt} = dum.savedata;
    end;
    
    aRT{it} = dat;
end;
clear dum dat;
%%
n = zeros(length(aRT),1);
parfor it = 1:length(aRT)
    
    n(it) = length(aRT{it});
    
end;

aRT2 = cell(sum(n),1);
idx = 1:length(aRT{1});
for it = 1:length(aRT)
    aRT2(idx) = aRT{it};
    if sign(length(aRT)-it)==1
        idx = idx(end)+1:idx(end)+length(aRT{it+1});
    end;
end;
clear aRT;
%%
idx = 1:size(RTdata{1},1);
rt = zeros(sum(n),1);
for it = 1:length(RTdata)
    
    rt(idx) = RTdata{it}(:,2);
    if sign(length(RTdata)-it)==1
        idx = idx(end)+1:idx(end)+size(RTdata{it+1},1);
    end;
end;
%%
del_idx = find(sign(rt-300)==-1);
aRT2(del_idx) = [];
rt(del_idx) = [];
%%
clear del_idx;
dat = aRT2; clear aRT*;
Y = rt; clear rt;
%% load the optimization data
path2files = '/bcbl/home/home_a-f/froux/BilingualNamingExp/greedy_optim/';
load([path2files,'greedy_optim_21-Aug-2015.mat']);

x = squeeze(median(optim_data.hist_t(:,end,:),1));
tresh = cell(6,1);
for it = 1:length(x)
    tresh{it} =x(it);
end;
on = zeros(length(dat),1);
parfor it = 1:length(dat)
    [on(it)] = detect_speech_on_and_offset(dat{it},tresh);
end;
%%
path2files = '/bcbl/home/public/YAZ/Expe_NP/NP/Mono_copy/audio_files/';
x = dir([path2files,'*.WAV']);

k =0;
ID2 = cell(sum(n),1);
for it = 1:length(files)
    for jt = 1:length(RTdata{it})
        
        chckf = dir([path2files,files(it).name(1:regexp(files(it).name,'.azk')-1),'_*',num2str(RTdata{it}(jt,1)),'.WAV']);
        if ~isempty(chckf)
            k = k+1;
            ID2{k} = chckf.name;
        end;
    end;
end;
%%
path2files = '/bcbl/home/public/YAZ/Expe_NP/NP/Mono_copy/audio_files/';

on2 = zeros(length(ID2),1);
parfor it = 1:length(ID2)
   
    
    wav_file = dir([path2files,ID2{it}]);
    wav_file.name = [path2files,wav_file.name];
    [on2(it)] = finalbatch_delay(wav_file);
    
end;
%%
X = 1:optim_data.Niter;
X = [1 X+1];

mlt = zeros(size(optim_data.hist_e,1),size(optim_data.hist_e,2)+1);
mlt(:,1) = optim_data.ei_training;
mlt(:,2:end) = optim_data.hist_e;

M = squeeze(mean(mlt,1));
SD = squeeze(std(mlt,1));


mlt = zeros(size(optim_data.test_e,1),size(optim_data.test_e,2)+1);
mlt(:,1) = optim_data.ei_test;
mlt(:,2:end) = optim_data.test_e;

M2 = squeeze(mean(mlt,1));
SD2 = squeeze(std(mlt,1));


figure;
subplot(321);
a = gca;
hold on;

shadedErrorBar(X,M,SD,{'Color',[.5 .5 .5]});

h = [];
h(1) = plot(X,M,'k','LineWidth',3,'LineSmoothing','on');
h(2) = plot(X,M2,'r','LineWidth',3,'LineSmoothing','on');

legend(h,'Training','Test');
ylim([min(M)-15 max(M)]);
xlabel('Iterations','Fontsize',14);
ylabel('MLE [\sigma]','Fontsize',14);


subplot(322);
a = [a gca];
hold on;

M1 = optim_data.hist_e(:,end);
M2 = optim_data.test_e(:,end);

D = M1-M2;

[n,x] = hist(D,[-50:5:35]);
[w,f] = ksdensity(D,[-50:5:50]);

bar(x,n./sum(n),'FaceColor','k');
plot(f,w./sum(w),'r','LineWidth',3,'LineSmoothing','on');

xlabel('[ms]','Fontsize',14);
ylabel('[PDF]','Fontsize',14);

subplot(323);
a = [a gca];
hold on;
plot(X(1:end-1),squeeze(mean(log10(optim_data.hist_o),1)),'k','LineWidth',3,'LineSmoothing','on');
xlabel('Iterations','Fontsize',14);
ylabel('\Omega [log]','Fontsize',14);

subplot(324);
a = [a gca];
hold on;

M = squeeze(mean(optim_data.hist_t(:,end,:),1));
SD = squeeze(std(optim_data.hist_t(:,end,:),1));

xlim([0 8]);

shadedErrorBar(1:6,M,SD,{'Color',[.9 0 0]});
h = [];
h(2) = plot(1:6,M,'rs-','LineWidth',3,'LineSmoothing','on');
h(1) = plot(1:6,[optim_data.itresh{:}],'ks-','LineWidth',3,'LineSmoothing','on');
set(gca,'XTick',[1:6]);
set(gca,'XTickLabel',{'Amp.' 'WE' 'SC' 'AM' 'FM' 'HP'},'Fontsize',14);
xlabel('Features','Fontsize',14);
ylabel('Threshold value [a.u.]','Fontsize',14);
legend(h,'Init.','Optim.');

% subplot(425);
% a = [a gca];
% hold on;
% 
mY = Y;
Lm = [min(on) max(on) min(mY) max(mY)];
Lm = [min(Lm) max(Lm)];
Lm = [floor(min(Lm)/1000)*1000 round(max(Lm)/1000).*1000];
% 
% [b1,bint,r,rint,stats1] = regress(Y(:,1),[ones(size(Y(:,1))) Y(:,2)]);
% plot(Lm,Lm,'r','LineWidth',3,'LineSmoothing','on');
% plot(Y(:,1),Y(:,2),'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
% axis tight;
% xlim([Lm]);
% ylim([Lm]);
% xlabel('Manual [ms]','Fontsize',14);
% ylabel('Manual [ms]','Fontsize',14);
% 
% text(550,2700,['R² = ',num2str(round(stats1(1)*1000)/1000)],'Fontsize',14);
% text(550,2525,['\alpha = ',num2str(round(b1(1)*1000)/1000)],'Fontsize',14);
% text(550,2350,['\beta = ',num2str(round(b1(2)*1000)/1000)],'Fontsize',14);

subplot(326);
a = [a gca];
hold on;

[b2,bint,r,rint,stats2] = regress(mY,[ones(size(mY(:,1))) on]);
plot(Lm,Lm,'r','LineWidth',3,'LineSmoothing','on');
plot(on,mY,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
axis tight;
xlim([Lm]);
ylim([Lm]);
xlabel('Automatic scores [ms]','Fontsize',14);
ylabel('Median manual scores [ms]','Fontsize',14);


text(550,2700,['R² = ',num2str(round(stats2(1)*1000)/1000)],'Fontsize',14);
text(550,2525,['\alpha = ',num2str(round(b2(1)*1000)/1000)],'Fontsize',14);
text(550,2350,['\beta = ',num2str(round(b2(2)*1000)/1000)],'Fontsize',14);

% subplot(427);
% a = [a gca];
% hold on;
% 
% [b3,bint,r,rint,stats3] = regress(mY,[ones(size(on2)) on2]);
% plot(Lm,Lm,'r','LineWidth',3,'LineSmoothing','on');
% plot(on2,mY,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
% axis tight;
% xlim([Lm]);
% ylim([Lm]);
% xlabel('Automatic scores [ms]','Fontsize',14);
% ylabel('Median manual scores [ms]','Fontsize',14);
% 
% 
% text(550,2700,['R² = ',num2str(round(stats3(1)*1000)/1000)],'Fontsize',14);
% text(550,2525,['\alpha = ',num2str(round(b3(1)*1000)/1000)],'Fontsize',14);
% text(550,2350,['\beta = ',num2str(round(b3(2)*1000)/1000)],'Fontsize',14);
% 
% subplot(428);
% a = [a gca];
% hold on;
% 
% [b4,bint,r,rint,stats4] = regress(mY,[ones(size(on3)) on3]);
% plot(Lm,Lm,'r','LineWidth',3,'LineSmoothing','on');
% plot(on3,mY,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
% axis tight;
% xlim([Lm]);
% ylim([Lm]);
% xlabel('Automatic scores [ms]','Fontsize',14);
% ylabel('Median manual scores [ms]','Fontsize',14);
% 
% text(550,2700,['R² = ',num2str(round(stats4(1)*1000)/1000)],'Fontsize',14);
% text(550,2525,['\alpha = ',num2str(round(b4(1)*1000)/1000)],'Fontsize',14);
% text(550,2350,['\beta = ',num2str(round(b4(2)*1000)/1000)],'Fontsize',14);
% 
% set(a(5:8),'XTick',Lm);
% set(a(5:8),'YTick',Lm);

axis(a,'tight');
set(a(4),'XLim',[0 8]);
set(a([1 3]),'XLim',[0 500]);
set(a,'Fontsize',14);
set(a,'LineWidth',3);
set(gcf,'Color','w');
%%
matlabpool('close');
%%
% path2files = '/bcbl/home/home_a-f/froux/BilingualNamingExp/manual_scores/';
% files = dir([path2files,'*.xlsx']);
% 
% [data3] = zeros(length(files),176,7);
% for it = 1:length(files)
%     [xdata,txt,raw] = xlsread([path2files,files(it).name]);
%     if size(xdata,1) == 177
%         xdata(end,:) = [];
%     end;
%     
%     data3(it,:,:) = xdata;
% end;
% clear xdata;
% data3(:,:,1) = [];
% 
% [data3] = data3(:,:,[1 3 5]);
% 
% 
% del_idx = [8,44,50,55,57,139,175];
% data3(:,del_idx,:) = [];
% data3(:,end,:) = [];
% 
% del_idx2 = cell(size(data3,1),1);
% for it = 1:size(data3,1)
%     
%     [del_idx2{it},~] = find(isnan(squeeze(data3(it,:,:))));
%     
% end;
% del_idx2 = unique([del_idx2{:}]);
% data3(:,del_idx2,:) = [];
% 
% mY3 = squeeze(median(data3,1));
% mY3 = [mY3(:,1);mY3(:,2);mY3(:,3)];
% dum = squeeze(median(data3(:,:,1)))';
% if ~isequal(dum,mY3(1:length(dum)))
%     error('data does not match');
% end;
% mY3 = mY3.*1000;
% 
% path2files = '/bcbl/home/home_a-f/froux/BilingualNamingExp/automatic_scores/feature_data/';
% files1 = dir([path2files,'Itzal-bilingual_picture_naming1_sound_recording-*.mat']);
% files2 = dir([path2files,'01_Itzal_Day2-bilingual_picture_naming_sound_recording-*.mat']);
% files3 = dir([path2files,'01_itzal_day3-bilingual_picture_naming_sound_recording-*.mat']);
% 
% files1(del_idx) = [];
% files2(del_idx) = [];
% files3(del_idx) = [];
% files1(end) = [];
% files2(end) = [];
% files3(end) = [];
% files1(del_idx2) = [];
% files2(del_idx2) = [];
% files3(del_idx2) = [];
% files = [files1;files2;files3];
% clear files1 files2 files3;
% 
% dat2 = cell(length(files),1);
% parfor it = 1:length(files)
%     
%     dum = load([path2files,files(it).name],'savedata');
%     dat2{it} = dum.savedata;
% 
% end;

% on2 = zeros(length(dat2),1);
% parfor it = 1:length(dat2)
%     
%     [on2(it)] = detect_speech_on_and_offset(dat2{it},tresh);
%     
% end;
% 
% rtr = zeros(size(data3,1),size(data3,2)*3);
% for it = 1:size(data3,1)
%     dum = squeeze(data3(it,:,:));
%     rtr(it,:) = dum(:);
% end;
% rtr = rtr.*1000;
% 
% C = corr(rtr');
% C(C==1) = NaN;
% [i1,i2] = find(C == max(max(C)));
% 
% rix = unique([i1 i2]);
% 
% D = abs(rtr(rix(1),:)' - rtr(rix(2),:)');
% 
% sel_idx = find(D<=30);
% 
% on2 = on2(sel_idx);
% 
% mY3 = median([rtr(i1,:)' rtr(i2,:)'],2);
% mY3 = mY3(sel_idx);