%%
restoredefaultpath;
addpath(genpath('/bcbl/home/home_a-f/froux/chronset/'));
%%
if matlabpool('size')==0
    matlabpool('local');%matlabpool local;
end;
%%
mode = 'orig';
%% load optimized thresholds for BCN data
load('/bcbl/home/home_a-f/froux/chronset/thresholds/greedy_optim_NP_data_BCN_28-Sep-2015.mat');
%load /bcbl/home/home_a-f/froux/chronset/thresholds/greedy_optim_NP_data_SayWhen_20-Jun-2016.mat;
%%
dum = optim_data.test_e;
dum(optim_data.test_e==0) = NaN;
ind = find(dum== min(min(dum)));
[i1,i2] = ind2sub(size(optim_data.test_e),ind);
%%
i1 = min(unique(i1));
i2 = min(unique(i2));

[i1 i2]
%%
tresh_params = cell(size(optim_data.hist_t,3),1);
for it = 1:size(optim_data.hist_t,3)
    tresh_params(it) = {squeeze(optim_data.hist_t(max(i1),min(i2),it))};%get set of thresholds that generated minimal test error
end;
%% read in the manual scores
path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/';

[~,txt,raw] = xlsread([path2files,'manual_ratings_BCN.xls']);
%% get rid of header
txt(1,:) = [];
raw(1,:) = [];

[pc] = txt(:,1); % participant ID
[fc] = [raw(:,4)]; % file ID

% convert cell to double
fc2 = zeros(length(fc)-1,1);
parfor it = 1:length(fc)
    fc2(it) = fc{it};
end;

if length(pc) ~= length(fc2)
    error('wrong number of files');
end;

Y = cell2mat(raw(:,9:10));% manual scores for the 2 raters
%% generate ID information

%build file ID based on .wav file-name
path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/raw_data/';
ID = cell(length(fc2),1);
parfor it = 1:length(fc2)
    
    [chck] = dir([path2files,pc{it},'_',num2str(fc2(it)),'.WAV']);
    ID{it} = chck.name;
    
end;

%double check make sure there's no missed file-names
[sel_idx] = zeros(length(ID),1);
parfor it = 1:length(ID)
    file_name = dir([path2files,ID{it}]);
    if ~isempty(file_name)        
        sel_idx(it) = it;
    end;
end;
[sel_idx] = sel_idx(sel_idx ~= 0);

% check
if (length(sel_idx) ~= length(pc)) || (length(sel_idx) ~= length(fc2)) || (length(sel_idx) ~= length(Y))
    error('Wrong number of files');
end;

% build file-name info based on .mat extension
switch mode
    case 'orig'
        path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/feature_data/';
    case 'new'
        path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/feature_data/Jun2016/';
end;

ID2 = cell(length(fc2),1);
for it = 1:length(fc2)
    
    switch mode
        case 'orig'
            [chck] = dir([path2files,pc{it},'_',num2str(fc2(it)),'*.mat']);
        case 'new'
            [chck] = dir([path2files,pc{it},'_',num2str(fc2(it)),'_*-Jun-2016.mat']);
    end;
    
    ID2{it} = chck.name;
    
end;

%check for missing files
[sel_idx] = zeros(length(ID2),1);
parfor it = 1:length(ID2)
    file_name = dir([path2files,ID2{it}]);
    if ~isempty(file_name)        
        sel_idx(it) = it;
    end;
end;
[sel_idx] = sel_idx(sel_idx ~= 0);

% error checks
if (length(sel_idx) ~= length(pc)) || (length(sel_idx) ~= length(fc2)) || (length(sel_idx) ~= length(ID2))
    error('Wrong number of files');
end;

chck = zeros(length(ID2),1);
parfor it = 1:length(ID2)    
    switch mode
        case 'orig'
            chck(it) = strcmp(ID{it}(1:end-4),ID2{it}(1:end-4));
        case 'new'
            chck(it) = strcmp(ID{it}(1:end-4),ID2{it}(1:max(regexp(ID2{it},'_')-1)));
    end;
end;
if sum(chck) ~= length(ID2)
    error('number of files does not match');
end;
%% compute onsets using chronset for BCN data based on optimized thresholds
switch mode
    case 'orig'
        path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/feature_data/';
    case 'new'
        path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/feature_data/Jun2016/';
end;

dat = cell(length(ID),1);

switch mode
    case 'orig'       
        
        on1 = zeros(length(dat),1);
        parfor it = 1:length(ID)
            
            fprintf([num2str(it),'/',num2str(length(ID))]);
            chck = dir([path2files,ID2{it}]);
            dum = load([path2files,chck.name]);
            [on1(it)] =  detect_speech_on_and_offset_orig2(dum.savedata,[tresh_params' {0.035} {4}]);
            fprintf('\n');
            
        end;
        clear dum;
    case 'new'
        idx1 = tresh_params{9};
        idx2 = tresh_params{10};
        
        on1 = zeros(length(dat),1);
        parfor it = 1:length(ID)
            
            fprintf([num2str(it),'/',num2str(length(ID))]);
            chck = dir([path2files,ID2{it}]);
            dum = load([path2files,chck.name]);
            [on1(it)] =  detect_speech_on_and_offset_orig2(dum.savedata{idx1,idx2},tresh_params);
            fprintf('\n');
            
        end;
        clear dum;
end;
%% write Chronset and human rater onset times for BCN to txt file
sfp = '/bcbl/home/home_a-f/froux/chronset/data/BCN/';
fid = fopen([sfp,'Chronset_onsets_vs_manual_onsets_BCNdata.txt'],'w+');

for it = 1:length(ID2)
    switch mode
        case 'orig'
            id = {ID2{it}(1:regexp(ID2{it},'*.mat')-1)};
        case 'new'
            id = {ID2{it}(1:regexp(ID2{it},'_*-Jun-2016.mat')-1)};
    end;
    id = [id{:},'.WAV'];
    y = [on1(it) Y(it,1) Y(it,2)];
    fprintf(fid,'%s\t',id);
    fprintf(fid,'%d\t%d\t%d\n',y);
end;
%% compute onsets using griffith method for BCN data (EDP)
path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/raw_data/';

on2 = zeros(length(ID2),1);
parfor it = 1:length(ID2)
    
    wav_file = [];
    switch mode
        case 'orig'
            wav_file = dir([path2files,ID{it}(1:end-4),'.WAV']);
        case 'new'
            wav_file = dir([path2files,ID{it}(1:max(regexp(ID2{it},'_'))-1),'.WAV']);
    end;
    
    wav_file.name = [path2files,wav_file.name];
    [on2(it)] = finalbatch_delay(wav_file);
    
end;
%% read in check volca based onset detection of BCN data
path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/';
file_name = [path2files,'check-output.txt'];
file = dir(file_name);

fid = fopen([path2files,file.name],'r');

[data] = textscan(fid,'%s');
data = data{:}';
k = 0;
data2 = cell(length(data)/2,2);
for it = 1:2:length(data)
    k = k+1;
    data2{k,1} = [data{it}];
    data2{k,2} = str2double([data{it+1}]);
end;

idx = zeros(length(ID2),1);
for it = 1:length(ID2)
    switch mode
        case 'orig'
            if ~isempty(find(strcmp([data2(:,1)],[ID2{it}(1:end-4),'.WAV'])))
                idx(it) = find(strcmp([data2(:,1)],[ID2{it}(1:end-4),'.WAV']));
            end;
        case 'new'
            if ~isempty(find(strcmp([data2(:,1)],[ID2{it}(1:max(regexp(ID2{it}(1:end-4),'_'))-1),'.WAV'])))
                idx(it) = find(strcmp([data2(:,1)],[ID2{it}(1:max(regexp(ID2{it}(1:end-4),'_'))-1),'.WAV']));
            end;
    end;
end;
idx = idx(find(idx~=0));
if length(idx) ~= length(ID2)
    error('wrong number of files detected');
end;

data2 = data2(idx,:);

[on3] = [data2{:,2}]';
%%
%Y = cell2mat(raw(:,9:10));% manual scores for the 2 raters
%mY1 = mean(Y,2);
[u,s,v] = svd(Y);
mY1 = u(:,1)*s(1,1)*v(:,1)';
%mY1 = mean(mY1,2);
%%
d1 = mY1(:,2)-mY1(:,1);
length(find(sign(10-abs(d1))==1))/length(d1)

d2 = mean(mY1,2)-on1;
length(find(sign(10-abs(d2))==1))/length(d2)

d3 = mean(mY1,2)-on2;
length(find(sign(10-abs(d3))==1))/length(d3)

d4 = mean(mY1,2)-on3;
length(find(sign(10-abs(d4))==1))/length(d4)
%%
d1 = (on1-mean(mY1,2));
d2 = (on2-mean(mY1,2));
d3 = (on3-mean(mY1,2));

[n1,x1] = hist(d1,[-200:200]);
[n2,x2] = hist(d2,[-200:200]);
[n3,x3] = hist(d3,[-200:200]);

p1 = sum(n1(find(x1 >= -10 & x1 <=10)))/sum(n1)
p2 = sum(n2(find(x2 >= -10 & x2 <=10)))/sum(n2)
p3 = sum(n3(find(x3 >= -10 & x3 <=10)))/sum(n3)
%% prepare data for making abs-Diff Histograms for BCN data
bs = [5:5:25 50:25:150 151];
bs = sort([-bs bs]);

[n1,x1] = hist((on1-mean(mY1,2)),bs);
[n2,x2] = hist((on2-mean(mY1,2)),bs);
[n3,x3] = hist((on3-mean(mY1,2)),bs);

hist_dat.onBCN1 = [x1' n1'];
hist_dat.onBCN2 = [x2' n2'];
hist_dat.onBCN3 = [x3' n3'];
%%
X = Y;
MSb = zeros(1,size(X,1));
for kt = 1:size(X,1)
    MSb(kt) = sum((mean(X,2)-mean(X(kt,:),2)).^2)/(size(X,1)-1);
end;

MSw = sum(sum((X-repmat(mean(X,2),[1 size(X,2)])).^2,2),1)/(size(X,1)*(size(X,2)-1));

ICC = (MSb-MSw)/(MSb+(size(X,2)-1)*MSw);
d1 = [X(:,1)-X(:,2)];

X = mY1;
MSb = zeros(1,size(X,1));
for kt = 1:size(X,1)
    MSb(kt) = sum((mean(X,2)-mean(X(kt,:),2)).^2)/(size(X,1)-1);
end;

MSw = sum(sum((X-repmat(mean(X,2),[1 size(X,2)])).^2,2),1)/(size(X,1)*(size(X,2)-1));

ICC2 = (MSb-MSw)/(MSb+(size(X,2)-1)*MSw);
d2 = [X(:,1)-X(:,2)];

mY1 = mean(Y,2);
[b3,bint,r1,rint,stats3] = regress(mY1,[ones(size(on1)) on1]);
[b4,bint,r2,rint,stats4] = regress(mY1,[ones(size(on2)) on2]);
[b5,bint,r3,rint,stats5] = regress(mY1,[ones(size(data2(:,1))) on3]);
[b6,bint,r6,rint,stats6] = regress(Y(:,2),[ones(size(Y(:,1))) Y(:,1)]);

Lm = [min([on1;on2;[data2{:,2}]']) max([on1;on2;[data2{:,2}]'])];

figure;
subplot(421);
a = gca;
hold on;
plot([Lm(1) Lm(2)],[Lm(1) Lm(2)],'r','LineWidth',3);
plot(Y(:,1),Y(:,2),'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
text(500,2500,{['a:',num2str(round(b6(1)*100)/100)];['R^2:',num2str(round(stats6(1)*100)/100)]});

subplot(422);
a = [a gca];
hold on;
plot([Lm(1) Lm(2)],[Lm(1) Lm(2)],'r','LineWidth',3);
plot(on1,mY1,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
text(500,2500,{['a:',num2str(round(b3(1)*100)/100)];['R^2:',num2str(round(stats3(1)*100)/100)]});

subplot(423);
a = [a gca];
hold on;
plot([Lm(1) Lm(2)],[Lm(1) Lm(2)],'r','LineWidth',3);
plot(on2,mY1,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
text(500,2500,{['a:',num2str(round(b4(1)*100)/100)];['R^2:',num2str(round(stats4(1)*100)/100)]});

subplot(424);
a = [a gca];
hold on;
plot([Lm(1) Lm(2)],[Lm(1) Lm(2)],'r','LineWidth',3);
plot([data2{:,2}]',mY1,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
text(500,2500,{['a:',num2str(round(b5(1)*100)/100)];['R^2:',num2str(round(stats5(1)*100)/100)]});

axis(a,'tight');
set(gcf,'Color','w');
set(a,'Xlim',Lm);
set(a,'Ylim',Lm);
set(a,'XTick',Lm);
set(a,'YTick',Lm);
set(a,'LineWidth',3);

xlabel(a(1),'Manual [ms]');
xlabel(a(2),'Chronset [ms]');
xlabel(a(3),'EdP [ms]');
xlabel(a(4),'Checkvocal [ms]');

for it = 1:length(a)   
    ylabel(a(it),'Avg. Manual [ms]');    
end;
    
R = [r1 r2 r3 r6];
bx = min(min(R)):1:max(max(R));


[n1,x1] = hist(r1,bx);
[n2,x2] = hist(r2,bx);
[n3,x3] = hist(r3,bx);
[n6,x6] = hist(r6,bx);

n1 = n1./sum(n1);
n2 = n2./sum(n2);
n3 = n3./sum(n3);
n6 = n6./sum(n6);

g1 = cumsum(n1);
g2 = cumsum(n2);
g3 = cumsum(n3);
g6 = cumsum(n6);

[idx1] = find(g1 >= .10 & g1 <=.9);
[idx2] = find(g2 >= .10 & g2 <=.9);
[idx3] = find(g3 >= .10 & g3 <=.9);

subplot(4,2,[5]);
hold on;
area([x1(min(idx1)) x1(max(idx1))],[1 1],0);
h = findobj('Type','patch');
set(h,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75]);

h = [];
h(1) = plot(x1,g1,'b','LineWidth',2);
h(2) = plot(x2,g2,'r','LineWidth',2);
plot([x1(idx1(1)) x1(idx1(1))],[0 1],'--','Color',[.75 .75 .75]);
plot([x1(idx1(end)) x1(idx1(end))],[0 1],'--','Color',[.75 .75 .75]);
plot([x1(idx1(1)) x1(idx1(end))],[0 0],'--','Color',[.75 .75 .75]);
plot([x1(idx1(1)) x1(idx1(end))],[1 1],'--','Color',[.75 .75 .75]);

plot([x2(idx2(1)) x2(idx2(1))],[0 1],'--','Color',[0 0 0]);
plot([x2(idx2(end)) x2(idx2(end))],[0 1],'--','Color',[0 0 0]);
plot([x2(idx2(1)) x2(idx2(end))],[0 0],'--','Color',[0 0 0]);
plot([x2(idx2(1)) x2(idx2(end))],[1 1],'--','Color',[0 0 0]);

axis tight;
ylim([0 1.05]);
xlabel('Regression residuals [ms]');
ylabel('CDF');
set(gca,'XTick',round(sort([x1(1) x2(idx2(1)) 0 x2(idx2(end)) 700])*1)/1);
set(gca,'YTick',[0 1]);
set(gca,'LineWidth',3);
set(gca,'Xlim',[x1(1) 700]);

subplot(4,2,[6]);
hold on;
area([x1(min(idx1)) x1(max(idx1))],[1 1],0);
h = findobj('Type','patch');
set(h,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75]);

h = [];
h(1) = plot(x1,g1,'LineWidth',2);
h(2) = plot(x3,g3,'r','LineWidth',2);
plot([x1(idx1(1)) x1(idx1(1))],[0 1],'--','Color',[.75 .75 .75]);
plot([x1(idx1(end)) x1(idx1(end))],[0 1],'--','Color',[.75 .75 .75]);
plot([x1(idx1(1)) x1(idx1(end))],[0 0],'--','Color',[.75 .75 .75]);
plot([x1(idx1(1)) x1(idx1(end))],[1 1],'--','Color',[.75 .75 .75]);

plot([x3(idx3(1)) x3(idx3(1))],[0 1],'--','Color',[0 0 0]);
plot([x3(idx3(end)) x3(idx3(end))],[0 1],'--','Color',[0 0 0]);
plot([x3(idx3(1)) x3(idx3(end))],[0 0],'--','Color',[0 0 0]);
plot([x3(idx3(1)) x3(idx3(end))],[1 1],'--','Color',[0 0 0]);

axis tight;
ylim([0 1.05]);
xlabel('Regression residuals [ms]');
ylabel('CDF');
set(gca,'XTick',round(sort([x1(1) x3(idx3(1)) 0 x3(idx3(end)) 700])*1)/1);
set(gca,'YTick',[0 1]);
set(gca,'LineWidth',3);
set(gca,'Xlim',[x1(1) 700]);

subplot(4,2,7);
joint_distribution_plot_measurement_errors(r1,r2,R);
xlabel('Chronset residuals [ms]');
ylabel('EdP residuals [ms]');
set(gca,'XTick',round([x1(1) x1(idx1(1)) x1(idx1(end)) 700].*1)/1)
set(gca,'YTick',round([x2(1) x2(idx2(1)) x2(idx2(end)) 700].*1)/1)
ylim([-292 700]);
xlim([-292 700]);

subplot(4,2,8);
joint_distribution_plot_measurement_errors(r1,r3,R);
xlabel('Chronset residuals [ms]');
ylabel('Checkvocal residuals [ms]');
set(gca,'XTick',round([x1(1) x1(idx1(1)) x1(idx1(end)) 700].*1)/1)
set(gca,'YTick',round([x3(1) x3(idx3(1)) x3(idx3(end)) 700].*1)/1)
ylim([-292 700]);
xlim([-292 700]);
 

p1 = length(find(r1 >= -10 & r1 <=10))/length(r1)
p2 = length(find(r2 >= -10 & r2 <=10))/length(r2)
p3 = length(find(r3 >= -10 & r3 <=10))/length(r3)
p4 = length(find(r6 >= -10 & r6 <=10))/length(r6)

std(r6)
std(r1)
std(r2)
std(r3)
%%
SIds = unique(pc);
rsq = zeros(length(SIds),1);
n = zeros(length(SIds),1);

h = figure(2);
xw = 8.5;
yw = 17.2;
for it = 1:length(SIds)
    idx = find(strcmp(pc,SIds(it)));
    n(it) = length(idx);
    ms = mY1(idx);
    x = on1(idx);            
    [~,~,~,~,stats] = regress(ms,[ones(size(x)) x]);
    rsq(it) = stats(1);
end;
[~,sidx] = sort(rsq,'descend');
lm1 = zeros(length(sidx),2);
lm2 = zeros(length(sidx),2);
a = zeros(length(sidx),1);
for it = 1:length(sidx)
    
    idx = find(strcmp(pc,SIds(sidx(it))));
    n(it) = length(idx);
    ms = mY1(idx);
    x = on1(idx);
    lm = [max(on1) max(mY1)];
    
    subplot(5,3,it+1);
    a(it) = gca;
    hold on;
    plot([0 max(lm)],[0 max(lm)],'r-','LineWidth',3);

    plot(x,ms,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
    axis tight;
    text(200,lm(2)-100,['R^2:',num2str(round(rsq(sidx(it))*100)/100)]);
    lm1(it,:) = get(gca,'XLim');
    lm2(it,:) = get(gca,'YLim');
end;
set(a,'XTick',[0 max(max(lm1))]);
set(a,'YTick',[0 max(max(lm2))]);
set(a,'XTick',[]);
set(a,'YTick',[]);

set(a,'LineWidth',3);
set(gcf,'Color','w');
set(a,'LineWidth',3);
set(a,'XLim',[0 max(max(lm1))]);
set(a,'YLim',[0 max(max(lm2))]);
set(a,'XTick',[]);
set(a,'YTick',[]);

subplot(5,3,1);
set(gca,'LineWidth',3);
set(gca,'XLim',[0 max(max(lm1))]);
set(gca,'YLim',[0 max(max(lm2))]);
set(gca,'XTick',round([0 max(max(lm1))]));
set(gca,'YTick',round([0 max(max(lm2))]));

set(h,'PaperSize',[8.5 17.21],'PaperUnits','centimeters','PaperPosition',[0,0,xw,yw]);
set(h,'PaperPositionMode','auto','PaperOrientation','portrait');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAYWHEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read the Checkvocal data
fid = fopen('/bcbl/home/home_a-f/froux/chronset/data/SayWhen/CheckVocal_SayWhen.txt','r');
CVdat = textscan(fid,'%s');
CVdat = CVdat{:}';
CVdat2 = cell(length(CVdat)/10,10);
idx = 1:10;
for it = 1:length(CVdat)/10
    CVdat2(it,:) = CVdat(idx);
    idx = idx+10;
end;
CVdat = CVdat2;clear CVdat2;

[del_idx,~] = find(isnan([CVdat{:,[5 6 8 9 10]}]));
del_idx = unique(del_idx);
CVdat(del_idx,:) = [];
CVdat(1,:) = [];

% get the onset ratings 5-6 BCBL huamns, 7 VoiceKey, 8 SaWhen human, 9 SayWhen algo, Checkvocal
dum= CVdat(:,[5 6 8 9 10]);
dum = str2double(dum);
%% generate file ID based on Excel data
switch mode
    case 'orig'
        path2featfiles = '/bcbl/home/home_a-f/froux/chronset/data/SayWhen/feature_data/';
        
    case 'new'
        path2featfiles = '/bcbl/home/home_a-f/froux/chronset/data/SayWhen/feature_data/Jun2016/';
        
end;

%get info from Excel data
ID = {CVdat{:,3}};% file id
trl_n = str2double(CVdat(:,4));% trl number

if size(dum,1) ~= length(ID)% check for inconsistency
    error('number of files must match');
end;

% go into feature dir and find files that match Excel data
ID2 = cell(length(ID),1);
sel_idx = zeros(length(ID),1);
k = 0;
for jt = 1:length(ID)
    
    % check if file exists
    switch mode
        case 'orig'
            chck = dir([path2featfiles,CVdat{jt,3},'.',num2str(CVdat{jt,4}),'*.mat']);
            if ~isempty(chck)
                [id] = [CVdat{jt,3},'.',num2str(CVdat{jt,4}),'.mat'];
                chck = dir([path2featfiles,id]);
            else
                [id] = [];
                [chck] = [];
            end;
        case 'new'
            chck = dir([path2featfiles,CVdat{jt,3},'.',num2str(CVdat{jt,4}),'_*-Jun-2016.mat']);
            if ~isempty(chck)
                [id] = [CVdat{jt,3},'.',num2str(CVdat{jt,4}),'_',chck.date(1:11),'.mat'];
                chck = dir([path2featfiles,id]);
            else
                [id] = [];
                [chck] = [];
            end;
    end;
        
    
    % if file exists double check if it fits with Excel data
    if ~isempty(chck) && strcmp(id,chck.name)
        % save the feature file ID
        k = k+1;
        ID2(k) = ID(jt);
        sel_idx(k) = jt;
    else
    end;
end;
%%
% remove exceeding pre-allocated data
ID2(k+1:end) = [];
sel_idx(k+1:end) = [];
% only keep trl-info for existing matlab files
trl_n = trl_n(sel_idx);
% only keep Excel data for existing matlab files
dum = dum(sel_idx,:);
CVdat = CVdat(sel_idx,:);
%% make sure the files have been read in correctly
chck = strcmp(ID2,ID(sel_idx)');
if any(chck==0)
    error('error number of files does not match');
end;

for jt = 1:length(sel_idx)
    id = [CVdat{jt,3},'.',num2str(CVdat{jt,4}),'.mat'];
    id2 = [ID2{jt},'.',num2str(trl_n(jt)),'.mat'];
    
    if ~strcmp(id,id2)
        error('file assignment does not match');
    end;
end;
%%
del_idx = find(strcmp(ID2,'S29'));

sel_idx(del_idx) = [];
ID2(del_idx) = [];
trl_n(del_idx) = [];
dum(del_idx,:) = [];
CVdat(del_idx,:) = [];
% %% load optimized thresholds for BCN data
% %load('/bcbl/home/home_a-f/froux/chronset/thresholds/greedy_optim_SayWhen_18-Sep-2015.mat');
% load /bcbl/home/home_a-f/froux/chronset/greedy_optim_NP_data_BCN_16-Jun-2016.mat;
% [i1,i2] = find(optim_data.test_e == min(min(optim_data.test_e)));%minimized mle for test data
% tresh_params = cell(size(optim_data.hist_t,3),1);
% for it = 1:size(optim_data.hist_t,3)
%     tresh_params(it) = {squeeze(optim_data.hist_t(max(i1),min(i2),it))};%get set of thresholds that generated minimal test error
% end;
%% use chronset to estimate speech onset from feature data
on1 = zeros(length(ID2),1);

parfor jt = 1:length(ID2)
    
    chck = [];
    switch mode
        case 'orig'
            chck = dir([path2featfiles,ID2{jt},'.',num2str(trl_n(jt)),'.mat']);
            dat =load([path2featfiles,chck.name]);            
            
            % onsets lign up with Excel data
            [on1(jt)] = detect_speech_on_and_offset_orig2(dat.savedata,[tresh_params' {0.035} {4}]);
        case 'new'
            idx1 = tresh_params{9};
            idx2 = tresh_params{10};
            chck = dir([path2featfiles,ID2{jt},'.',num2str(trl_n(jt)),'_*-Jun-2016.mat']);
            dat =load([path2featfiles,chck.name]);
            dat = dat.savedata{idx1,idx2};
            
            % onsets lign up with Excel data
            [on1(jt)] = detect_speech_on_and_offset_orig2(dat,tresh_params);
    end;
    
end;
clear dum2;
%% use EPd to estimate speech onset from feature data
path2wavfiles = '/bcbl/home/home_a-f/froux/chronset/data/SayWhen/raw_data/';

on2 = zeros(length(ID2),1);
parfor jt = 1:length(ID2)
    
    [wav_file] =dir([path2wavfiles,ID2{jt},'.',num2str(trl_n(jt)),'.wav']);
    wav_file.name = [path2wavfiles,wav_file.name];
    % onsets lign up with Excel data
    [on2(jt)] = finalbatch_delay(wav_file);
end;
clear dum2;
%%
[on3] = str2double(CVdat(:,end));
%%
[on4] = dum(:,end-1);
%%
del_idx = [find(isnan(on1)) find(isnan(on2))];
on1(del_idx) = [];
on2(del_idx) = [];
on3(del_idx) = [];
on4(del_idx) = [];

dum(del_idx,:) = [];
ID2(del_idx,:) = [];
CVdat(del_idx,:) = [];
%% prepare data for making abs-Diff Histograms for SayWhen data
[u,s,v] = svd(dum(:,1:3));
mY1 = mean(u(:,1)*s(1,1)*v(:,1)',2);

%mY1 = mean(dum(:,1:3),2);

bs = [5:5:25 50:25:150 151];
bs = sort([-bs bs]);

[n1,x1] = hist((on1-mY1),bs);
[n2,x2] = hist((on2-mY1),bs);
[n3,x3] = hist((on3-mY1),bs);
[n4,x4] = hist((on4-mY1),bs);


hist_dat.onSW1 = [x1' n1'];
hist_dat.onSW2 = [x2' n2'];
hist_dat.onSW3 = [x3' n3'];
hist_dat.onSW4 = [x4' n4'];
%%
d1 = dum(:,1)-dum(:,2);
p1 = length(find(sign(10-abs(d1))==1))/length(d1)

d2 = dum(:,2)-dum(:,3);
p2 = length(find(sign(10-abs(d2))==1))/length(d2)

d3 = dum(:,1)-dum(:,3);
p3 = length(find(sign(10-abs(d3))==1))/length(d3)

d4 = mean(mY1,2)-on1;
p4 = length(find(sign(10-abs(d4))==1))/length(d4)

d5 = mean(mY1,2)-on2;
p5 = length(find(sign(10-abs(d5))==1))/length(d5)

d6 = mean(mY1,2)-on3;
p6 = length(find(sign(10-abs(d6))==1))/length(d6)

d7 = mean(mY1,2)-on4;
p7 = length(find(sign(10-abs(d7))==1))/length(d7)
%%
[b1,~,r1,~,stats1]= regress(dum(:,1),[ones(size(dum(:,2))) dum(:,2)]);
ypred1 = b1(1) +b1(2).*dum(:,2);

[b2,~,r2,~,stats2]= regress(dum(:,3),[ones(size(mean(dum(:,1:2),2))) mean(dum(:,1:2),2)]);
ypred2 = b2(1) +b2(2).*mean(dum(:,1:2),2);

Y = mean(dum(:,1:3),2);
[b3,~,r3,~,stats3]= regress(Y,[ones(size(on1)) on1]);%CHRON
ypred3 = b3(1)+b3(2).*on1;

[b4,~,r4,~,stats4]= regress(Y,[ones(size(on2)) on2]);%EDP
ypred4 = b4(1)+b4(2).*on2;

[b5,~,r5,~,stats5]= regress(Y,[ones(size(on3)) on3]);%CHECKV
ypred5 = b5(1)+b5(2).*on3;

on4= dum(:,4);
Y2 = mean(dum(:,3),2);
[b6,~,r6,~,stats6]= regress(Y2,[ones(size(on4)) on4]);%SW
ypred6 = b6(1)+b6(2).*on4;

lm = [max(max(dum)) max(on1) max(on2) max(on3) max(on4)];
lm = [0 max(lm)];
%%
figure;
subplot(431);
a = gca;
hold on;
plot([lm],[lm],'r','LineWidth',3);
%plot(dum(:,2),ypred1,'k','LineWidth',3);
plot(dum(:,2),dum(:,1),'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
text(1000,4500,{['b:',num2str(round(b1(2)*100)/100)];['a:',num2str(round(b1(1)*100)/100)];['R^2:',num2str(round(stats1(1)*100)/100)]});
xlabel('Manual [s]');
ylabel('Manual [s]');

subplot(432);
a = [a gca];
hold on;
plot([lm],[lm],'r','LineWidth',3);
%plot(mean(dum(:,1:2),2),ypred2,'k','LineWidth',3);
plot(mean(dum(:,1:2),2),dum(:,3),'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
text(1000,4500,{['b:',num2str(round(b2(2)*100)/100)];['a:',num2str(round(b2(1)*100)/100)];['R^2:',num2str(round(stats2(1)*100)/100)]});
xlabel('Manual [s]');
ylabel('Manual [s]');

subplot(433);
a = [a gca];
hold on;
plot([lm],[lm],'r','LineWidth',3);
%plot(on1,ypred3,'k','LineWidth',3);
plot(on1,Y,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
axis tight;
text(1000,4500,{['b:',num2str(round(b3(2)*100)/100)];['a:',num2str(round(b3(1)*100)/100)];['R^2:',num2str(round(stats3(1)*100)/100)]});
xlabel('Chronset [s]');
ylabel('Manual [s]');

subplot(434);
a = [a gca];
hold on;
plot([lm],[lm],'r','LineWidth',3);
%plot(on4,ypred6,'k','LineWidth',3);
plot(on4,Y2,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
axis tight;
text(1000,4500,{['b:',num2str(round(b6(2)*100)/100)];['a:',num2str(round(b6(1)*100)/100)];['R^2:',num2str(round(stats6(1)*100)/100)]});
xlabel('SayWhen [s]');
ylabel('Manual [s]');

subplot(435);
a = [a gca];
hold on;
plot([lm],[lm],'r','LineWidth',3);
%plot(on2,ypred4,'k','LineWidth',3);
plot(on2,Y,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
axis tight;
text(1000,4500,{['b:',num2str(round(b4(2)*100)/100)];['a:',num2str(round(b4(1)*100)/100)];['R^2:',num2str(round(stats4(1)*100)/100)]});
xlabel('Epd [s]');
ylabel('Manual [s]');

subplot(436);
a = [a gca];
hold on;
plot([lm],[lm],'r','LineWidth',3);
%plot(on3,ypred5,'k','LineWidth',3);
plot(on3,Y,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
axis tight;
text(1000,4500,{['b:',num2str(round(b5(2)*100)/100)];['a:',num2str(round(b5(1)*100)/100)];['R^2:',num2str(round(stats5(1)*100)/100)]});
xlabel('CheckVocal [s]');
ylabel('Manual [s]');

set(a,'LineWidth',3);
set(a,'XLim',[lm]);
set(a,'YLim',[lm]);
set(a,'XTick',[lm]);
set(a,'YTick',[lm]);
set(a,'XTickLabel',round([lm]./1000));
set(a,'YTickLabel',round([lm]./1000));

R = [r1 r2 r3 r4 r5 r6];

bx = round(round(min(min(R)))/10)*10:10:round(round(max(max(R)))/10)*10;
%bx = linspace(bx(1),bx(end),length(bx));

[n1,x1] = hist(r3,bx);
[n2,x2] = hist(r6,bx);
[n3,x3] = hist(r4,bx);
[n4,x4] = hist(r5,bx);

n1 = n1./sum(n1);
n2 = n2./sum(n2);
n3 = n3./sum(n3);
n4 = n4./sum(n4);

g1 = cumsum(n1);
g2 = cumsum(n2);
g3 = cumsum(n3);
g4 = cumsum(n4);

[idx1] = find(g1 >= 0.1 & g1 <= 0.9);
[idx2] = find(g2 >= 0.1 & g2 <= 0.9);
[idx3] = find(g3 >= 0.1 & g3 <= 0.9);
[idx4] = find(g4 >= 0.1 & g4 <= 0.9);

p1 = length(find(r3 >= -10 & r3 <=10))/length(r3)
p2 = length(find(r6 >= -10 & r6 <=10))/length(r6)
p3 = length(find(r4 >= -10 & r4 <=10))/length(r4)
p4 = length(find(r5 >= -10 & r5 <=10))/length(r5)

p5 = length(find(r1 >= -10 & r1 <=10))/length(r1)
p6 = length(find(r2 >= -10 & r2 <=10))/length(r2)

std(r3)
std(r6)
std(r4)
std(r5)



subplot(437);
a = gca;
hold on;
area([x1(min(idx1)) x1(max(idx1))],[1 1],0);
h = findobj('Type','patch');
set(h,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75]);

plot([x2(min(idx2)) x2(min(idx2))],[0 1],'k--');
plot([x2(max(idx2)) x2(max(idx2))],[0 1],'k--');
plot([x2(min(idx2)) x2(max(idx2))],[1 1],'k--');
plot([x2(min(idx2)) x2(max(idx2))],[0 0],'k--');

h = [];
h(2) = plot(x1,sort(g1)','b','LineWidth',3);
h(1) = plot(x2,sort(g2)','r','LineWidth',3);
axis tight;
ylim([0 1.05]);
ylabel('CDF');
xlabel('Residuals [ms]');
set(gca,'YTick',[0 1]);
set(gca,'LineWidth',3);
set(gca,'XTick',[x2(min(idx2)) x2(max(idx2)) 1.5e3]);

subplot(438);
a = [a gca];
hold on;
area([x1(min(idx1)) x1(max(idx1))],[1 1],0);
h = findobj('Type','patch');
set(h,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75]);

plot([x2(min(idx3)) x2(min(idx3))],[0 1],'k--');
plot([x2(max(idx3)) x2(max(idx3))],[0 1],'k--');
plot([x2(min(idx3)) x2(max(idx3))],[1 1],'k--');
plot([x2(min(idx3)) x2(max(idx3))],[0 0],'k--');

h = [];
h(2) = plot(x1,sort(g1)','b','LineWidth',3);
h(1) = plot(x2,sort(g3)','r','LineWidth',3);
axis tight;
ylim([0 1.05]);
ylabel('CDF');
xlabel('Residuals [ms]');
set(gca,'YTick',[0 1]);
set(gca,'LineWidth',3);
set(gca,'XTick',[x3(min(idx3)) x3(max(idx3)) 1.5e3]);

subplot(439);
a = [a gca];
hold on;
area([x1(min(idx1)) x1(max(idx1))],[1 1],0);
h = findobj('Type','patch');
set(h,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75]);

plot([x2(min(idx4)) x2(min(idx4))],[0 1],'k--');
plot([x2(max(idx4)) x2(max(idx4))],[0 1],'k--');
plot([x2(min(idx4)) x2(max(idx4))],[1 1],'k--');
plot([x2(min(idx4)) x2(max(idx2))],[0 0],'k--');

h = [];
h(2) = plot(x1,sort(g1)','b','LineWidth',3);
h(1) = plot(x4,sort(g4)','r','LineWidth',3);
axis tight;
ylim([0 1.05]);
ylabel('CDF');
xlabel('Residuals [ms]');
set(gca,'YTick',[0 1]);
set(gca,'LineWidth',3);
set(gca,'XTick',[x4(min(idx4)) x4(max(idx4)) 1.5e3]);

set(a,'XLim',[-750 1.5e3]);

subplot(4,3,10);
a = [gca];
joint_distribution_plot_measurement_errors(r3,r6,R);
xlabel('Chronset [ms]');
ylabel('SayWhen [ms]');
set(gca,'XTick',round([x1(idx1(1)) x1(idx1(end))].*1)/1)
set(gca,'YTick',round([x2(idx2(1)) x2(idx2(end))].*1)/1)
colorbar off;

subplot(4,3,11);
a = [a gca];
joint_distribution_plot_measurement_errors(r3,r4,R);
xlabel('Chronset [ms]');
ylabel('Epd [ms]');
set(gca,'XTick',round([x1(idx1(1)) x1(idx1(end))].*1)/1)
set(gca,'YTick',round([x3(idx3(1)) x3(idx3(end))].*1)/1)
colorbar off;

subplot(4,3,12);
a = [a gca];
joint_distribution_plot_measurement_errors(r3,r5,R);
xlabel('Chronset [ms]');
ylabel('Checkvocal [ms]');
set(gca,'XTick',round([x1(idx1(1)) x1(idx1(end))].*1)/1)
set(gca,'YTick',round([x4(idx4(1)) x4(idx4(end))].*1)/1)
colorbar off;

set(a,'YLim',[-1000 2500]);
set(a,'XLim',[-1000 2500]);

set(gcf,'Color','w');
%%
d1 = (on1-mY1);
d2 = (on4-mY1);
d5 = abs(d1);
d6 = abs(d2);

bw = linspace(-200,200,1000);

[n1,x1] = hist(r3,bw);
[n2,x2] = hist(r6,bw);
[n3,x3] = hist(d1,bw);
[n4,x4] = hist(d2,bw);

bw = linspace(0,200,500);

[n5,x5] = hist(d5,bw);
[n6,x6] = hist(d6,bw);

p1 = sum(n1(find(x1 >= -10 & x1 <=10)))/sum(n1)
p2 = sum(n2(find(x2 >= -10 & x2 <=10)))/sum(n2)
p3 = sum(n3(find(x3 >= -10 & x4 <=10)))/sum(n3)
p4 = sum(n4(find(x4 >= -10 & x4 <=10)))/sum(n4)




figure;
subplot(321);
hold on;
plot(x1,n1,'b');
plot(x2,n2,'r');
axis tight;
legend('Chronset','SayWhen');
title('Regression residuals');
subplot(323);
hold on;
plot(x3,n3,'b');
plot(x4,n4,'r');
axis tight;
title('Raw differences');
subplot(325);
hold on;
plot(x5,n5,'b');
plot(x6,n6,'r');
axis tight;
title('Abs. differences');

subplot(322);
hold on;
plot(x1,cumsum(n1)./sum(n1),'b');
plot(x2,cumsum(n2)./sum(n2),'r');
axis tight;
title('Regression residuals');
subplot(324);
hold on;
plot(x3,cumsum(n3)./sum(n3),'b');
plot(x4,cumsum(n4)./sum(n4),'r');
axis tight;
title('Raw differences');
subplot(326);
hold on;
plot(x5,cumsum(n5)./sum(n5),'b');
plot(x6,cumsum(n6)./sum(n6),'r');
axis tight;
title('Abs. differences');

set(gcf,'Color','w');
%%
d1 = (on1-mY1);
d2 = (on4-mY1);

d3 = abs(on1-mY1);
d4 = abs(on4-mY1);

bs = [5 10 15 25 50 100 150 250 500 750];
bs2 = [0:25:750];
bs3 = [0:50:750];
bs4 = [0:75:750];


[n1,x1] = compute_pct_trial(d3,bs);
[n2,x2] = compute_pct_trial(d4,bs);

[n3,x3] = compute_pct_trial(d3,bs2);
[n4,x4] = compute_pct_trial(d4,bs2);

[n5,x5] = compute_pct_trial(d3,bs3);
[n6,x6] = compute_pct_trial(d4,bs3);

[n7,x7] = compute_pct_trial(d3,bs4);
[n8,x8] = compute_pct_trial(d4,bs4);


figure;

subplot(4,1,1);
a = gca;
hold on;
bar(1:2:2*length(n1),n1,.35);
bar(2:2:2*length(n2),n2,.35,'r');
set(gca,'XTick',1.5:2:2*length(n1));
set(gca,'XTickLabel',bs);

subplot(4,1,2);
a = [a gca];
hold on;
bar(1:2:2*length(n3),n3,.35);
bar(2:2:2*length(n4),n4,.35,'r');
set(gca,'XTick',1.5:2:2*length(n3));
set(gca,'XTickLabel',bs2);

subplot(4,1,3);
a = [a gca;];
hold on;
bar(1:2:2*length(n5),n5,.35);
bar(2:2:2*length(n6),n6,.35,'r');
set(gca,'XTick',1.5:2:2*length(n5));
set(gca,'XTickLabel',bs3);

subplot(4,1,4);
a = [a gca];
hold on;
bar(1:2:2*length(n7),n7,.35);
bar(2:2:2*length(n8),n8,.35,'r');
set(gca,'XTick',1.5:2:2*length(n7));
set(gca,'XTickLabel',bs4);

for it = 1:length(a)
    xlabel(a(it),'Difference score [ms]');
    ylabel(a(it),'Percent trials [%]');
end;

set(gcf,'Color','w');
%%
figure;
SIds = unique(ID2);
rsq = zeros(length(SIds),1);
n = zeros(length(SIds),1);

lm = [max(max(dum)) max(on1) max(on2)];
lm = [0 max(lm)];
pc = cell(length(SIds),1);
for it = 1:length(SIds)
    
    idx = find(strcmp(ID2,SIds(it)));
    if unique(strcmp(SIds(it),CVdat(idx,3))) ==1
        n(it) = length(idx);
        pc(it) = unique(ID2(idx));
        ms = mean(dum(idx,1:2),2);
        x = on1(idx);
        
        [~,~,~,~,stats] = regress(ms,[ones(size(x)) x]);
        rsq(it) = stats(1);
    else
        error('wrong file assignment');
    end;
end;

[~,s_idx] = sort(rsq,'descend');
subplot(8,4,1);
axis tight;
set(gca,'XLim',round([100 max(lm)]));
set(gca,'YLim',round([100 max(lm)]));

set(gca,'XTick',round([100 max(lm)]));
set(gca,'YTick',round([100 max(lm)]));
set(gca,'LineWidth',3);

a = zeros(length(s_idx),1);
for it = 1:length(s_idx)       
    
    idx = find(strcmp(ID2,SIds(s_idx(it))));
    
    ms = mean(dum(idx,1:3),2);
    x = on1(idx);
   
    subplot(7,4,it+1);
    a(it) = gca;
    hold on;
    plot([0 max([max(on1) max(mean(dum(:,1:3),2))])],[0 max([max(on1) max(mean(dum(:,1:3),2))])],'r-','LineWidth',3);
    
    plot(x,ms,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
    axis tight;
    text(500,max([max(x) max(ms)])-100,num2str(round(rsq(s_idx(it))*100)/100));

    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    
end;
set(a,'Xlim',[0 max([max(on1) max(mean(dum(:,1:3),2))])]);
set(a,'YLim',[0 max([max(on1) max(mean(dum(:,1:3),2))])]);


set(gcf,'Color','w');
set(a,'LineWidth',3);
%%
visualize_results_MC_simulation;
set(gca,'YTick',[0:.1:1]);
ylim([-.01 .6]);
%%
[lds] = ds_comparison(on1,on4);

ds.ds1 = on1;% Chronset
ds.ds2 = on4;% SW
ds.lds = lds;% Chron-SW
ds.g = mY1;%manual

figure;
visualizeation_ds_comparions(ds);

return;
% %%
for it = 1:gcf
    set(it,'PaperPositionMode','auto');
end;
%%
savepath = '/bcbl/home/home_a-f/froux/chronset/figures/EPS/';
for it = 1:gcf
    print(it,'-depsc2','-adobecset','-paint','-r300',[savepath,'Chronset_figure_',num2str(10),'.eps']);
end;
% %%
% savepath = '/bcbl/home/home_a-f/froux/chronset/figures/EPS/';
% for it = 1:gcf
%     print(it,'-depsc2','-adobecset','-paint','-r300',[savepath,'Chronset_colorbar.eps']);
% end;   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    







