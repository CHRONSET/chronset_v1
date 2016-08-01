%%
restoredefaultpath;
addpath(genpath('/bcbl/home/home_a-f/froux/chronset/'));
%%
p2df = '/bcbl/home/home_a-f/froux/chronset/data/SayWhen/';
fn = 'file_and_phon.txt';

[pc] = read_phono_coding_info(p2df,fn);
%% read the SayWhen data
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
path2featfiles = '/bcbl/home/home_a-f/froux/chronset/data/SayWhen/feature_data/';

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
    [id] = [CVdat{jt,3},'.',num2str(CVdat{jt,4}),'.mat'];
    chck = dir([path2featfiles,id]);
    
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
%% load optimized thresholds for BCN data
%load('/bcbl/home/home_a-f/froux/chronset/thresholds/greedy_optim_SayWhen_18-Sep-2015.mat');
load /bcbl/home/home_a-f/froux/chronset/thresholds/greedy_optim_NP_data_BCN_28-Sep-2015.mat;
[i1,i2] = find(optim_data.test_e == min(min(optim_data.test_e)));%minimized mle for test data
tresh_params = cell(size(optim_data.hist_t,3),1);
for it = 1:size(optim_data.hist_t,3)
    tresh_params(it) = {squeeze(optim_data.hist_t(max(i1),min(i2),it))};%get set of thresholds that generated minimal test error
end;
%% open matlabpool
if matlabpool('size')==0
    matlabpool(32);
end;
%% use chronset to estimate speech onset from feature data
on1 = zeros(length(ID2),1);
%inf_f = zeros(length(ID2),6);

on_c = zeros(length(ID2),6);

parfor jt = 1:length(ID2)
    
    fprintf([num2str(jt),'/',num2str(length(ID2))]);
    dat =load([path2featfiles,ID2{jt},'.',num2str(trl_n(jt)),'.mat']);
    dat = dat.savedata;    
    
    % onsets lign up with Excel data
    [on1(jt),~,feat_dat] = detect_speech_on_and_offset_orig2(dat,tresh_params);
    
    %inf_f(jt,:) = feat_dat.finf.inf_f;
    
    x = zeros(6,length(feat_dat.finf.x_tresh{1}));
    for kt = 1:length(feat_dat.finf.x_tresh)
        x(kt,:) = feat_dat.finf.x_tresh{kt};
    end;
    
    [idx] = min(find(sign(sum(x,1)-3)==1));
        
    y = NaN(length(feat_dat.finf.x_tresh),1);
    
    if ~isempty(idx)
        for kt = 1:length(feat_dat.finf.x_tresh)
            %idx            
            y(kt) = feat_dat.finf.x_tresh{kt}(idx);           
        end;
    end;
    %y
    on_c(jt,:) = y;
    
    fprintf('\n');
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
%% close matlabpool
matlabpool close;
%% Checkvocal
[on3] = dum(:,end);%str2double(CVdat(:,end-1));
%% SayWhen
[on4] = dum(:,end-1);
%%
del_idx = [find(isnan(on1)) find(isnan(on2))];
on1(del_idx) = [];
on2(del_idx) = [];
on3(del_idx) = [];
on4(del_idx) = [];
dum(del_idx,:) = [];
data(del_idx,:) = [];
CVdat(del_idx,:) = [];
%%
id1 = CVdat(:,2);
for it = 1:length(id1)
    id1{it} = id1{it}(1:end-4);
end;
id2 = pc(:,1);
%%
CVdat(1,end+1) = {'0'};

x = zeros(length(id2),1);
k = 0;
for it = 1:length(id2)
    
    %lookup
    sel_idx = find(strcmp(id1(:),id2(it)));
    
    % if there is an entry
    if ~isempty(sel_idx)        
        CVdat(sel_idx,end) = pc(it,2);%x(k)                        
    end;
    
end;
o = zeros(size(CVdat,1),1);
for it = 1:size(CVdat,1);
    o(it) = isempty(CVdat{it,end});
end;
CVdat(find(o==1),:) = [];
%%
phID = unique(CVdat(:,end));
SA1 = cell(length(phID),1);
SA2 = cell(length(phID),1);
SA3 = cell(length(phID),1);
SA4 = cell(length(phID),1);

n = zeros(length(phID),1);
for it = 1:length(phID)
    
    sel_idx = find(strcmp(CVdat(:,end),phID(it)));
    
    n(it) = length(sel_idx);
    
    [u,s,v] = svd(dum(sel_idx,1:3));
    rec =  u(:,1)*s(1,1)*v(:,1)';
    rec = mean(rec,2);
    
    SA1{it} = [on1(sel_idx), rec];%Chronset (1)
    SA2{it} = [on2(sel_idx), rec];%Edp (2)
    SA3{it} = [on3(sel_idx), rec];%CheckVocal (3)
    SA4{it} = [on4(sel_idx), rec];%SayWhen (4)

end;
%%
n2 = cell(length(phID),1);

for it = 1:length(phID)
    x = (on_c(find(strcmp(phID(it),CVdat(:,end))),:));
    n2{it} = sum(x,1)./size(x,1);
end;
%%
a = zeros(length(n2),1);
figure;
k= 0;
for it = 1:length(n2)
    k = k+1;
    subplot(ceil(length(n2)/5),5,k);
    a(it) = gca;
    bar(1:6,n2{it});
    axis tight;
    xlim([0 7]);
    ylim([0 1]);
    title(phID(it));
end;
set(gcf,'Color','w');
for it = 1:length(a)
    xlabel(a(it),'Features');
    ylabel(a(it),'[%]');
end;
%%
figure;bar(1:length(phID),n);set(gca,'XTick',1:length(phID));set(gca,'XTickLabel',phID);
%%
% this section of the code does a hard-coded outlier screening and removes
% all onsets below 250 ms 
[params1] = compute_slope_parameters4onsets(SA1,phID,n);
[params2] = compute_slope_parameters4onsets(SA2,phID,n);
[params3] = compute_slope_parameters4onsets(SA3,phID,n);
[params4] = compute_slope_parameters4onsets(SA4,phID,n);

sel_idx1 = find(sign(params1.n-20)==1);
sel_idx2 = find(sign(params2.n-20)==1);
sel_idx3 = find(sign(params3.n-20)==1);
sel_idx4 = find(sign(params4.n-20)==1);

[params1] = sort_params(params1,sel_idx1);
[params2] = sort_params(params2,sel_idx2);
[params3] = sort_params(params3,sel_idx3);
[params4] = sort_params(params4,sel_idx4);
%%
lattice_plot_1(params1);
%lattice_plot_1(params2);
%lattice_plot_1(params3);
lattice_plot_1(params4);
%%
for it = 1:2%
    set(it,'PaperPositionMode','auto');
end;
%%
savepath = '/bcbl/home/home_a-f/froux/chronset/figures/EPS/';
for it = 1:2
    print(it,'-depsc2','-adobecset','-paint','-r300',[savepath,'Chronset_lattice_plot_phono_coding',num2str(it),'.eps']);
end;
%%
figure;
subplot(221);
a = gca;
cross_hair_plot(params3);
title('CheckVocal');
subplot(222);
a = [a gca];
cross_hair_plot(params2);
title('EdP');

xl = get(a,'XLim');
yl = get(a,'YLim');

set(a,'XLim',[min([xl{:}]) max([xl{:}])]);
set(a,'YLim',[min([yl{:}]) max([yl{:}])]);

subplot(223);
a = [gca];
cross_hair_plot(params4);
title('SayWhen');
subplot(224);
a = [a gca];
cross_hair_plot(params1);
title('Chronset');

xl = get(a,'XLim');
yl = get(a,'YLim');

set(a,'XLim',[min([xl{:}]) max([xl{:}])]);
set(a,'YLim',[min([yl{:}]) max([yl{:}])]);
set(gcf,'Color','w');
%%
figure;
plot(sort([params1.rsq-params4.rsq]).*1,'bo-','MarkerFaceColor','b');
axis tight;
set(gca,'XTick',1:length(params1.phID));
set(gca,'XTickLabel',params1.phID(:));

xlabel('Phonetic code');
ylabel('R² difference');
set(gcf,'Color','w');
set(gca,'Color','w');

wrs1 = sum(params1.rsq.*params1.n)/sum(params1.n)
wrs4 = sum(params4.rsq.*params4.n)/sum(params4.n)