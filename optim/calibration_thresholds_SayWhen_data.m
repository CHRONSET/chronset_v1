%%
restoredefaultpath;
addpath(genpath('~/froux/chronset/'));
%%
f = 1;
if f ==1
    if matlabpool('size')==0
        matlabpool local;%130;
    end;
    % number of partitions
    Ntrl = matlabpool('size');
elseif f ==0
    if matlabpool('size')==0
        matlabpool local;
    end;
    % number of partitions
    Ntrl = 1;
end;
%%
[data,txt,raw] = xlsread('~/froux/chronset/data/SayWhen/manual_ratings_SayWhen.xls');
raw(1,:) = [];
txt(1,:) = [];
% %%
% data([210,218,221,370,397,429,433,470,619,1320,1783,1956,2010,3039,3141,3193,3664,4191,4290,4840],:)=[];
% txt([210,218,221,370,397,429,433,470,619,1320,1783,1956,2010,3039,3141,3193,3664,4191,4290,4840],:)=[];
% raw([210,218,221,370,397,429,433,470,619,1320,1783,1956,2010,3039,3141,3193,3664,4191,4290,4840],:)=[];
%%
[del_idx] = find(strcmp(raw(:,2),'..wav'));
raw(del_idx,:) = [];
txt(del_idx,:) = [];
data(del_idx,:) = [];
%%
[del_idx] = [find(isnan(data(:,5)));find(isnan(data(:,6)));find(isnan(data(:,7)));find(isnan(data(:,8)));find(isnan(data(:,9)))];
raw(del_idx,:) = [];
txt(del_idx,:) = [];
data(del_idx,:) = [];
%%
dum= data(:,[5:9]);
[U,S,V] = svd(dum(:,[1,2,4]));
pred = U(:,1)*S(1,1)*V(:,1)';
%%
path2featfiles = '~/froux/chronset/data/SayWhen/feature_data/';

ID = {txt{:,3}};
trl_n = [raw{:,4}];

if size(pred,1) ~= length(ID)
    error('number of files must match');
end;

ID2 = cell(length(ID),1);
sel_idx = zeros(length(ID),1);
k = 0;
for jt = 1:length(ID)
    
    chck = dir([path2featfiles,ID{jt},'.',num2str(trl_n(jt)),'.mat']);
    [id] = [raw{jt,3},'.',num2str(raw{jt,4}),'.mat'];    
    
    if ~isempty(chck) && strcmp(id,chck.name)
        k = k+1;
        ID2(k) = ID(jt);
        sel_idx(k) = jt;
    end;
end;
ID2(k+1:end) = [];
sel_idx(k+1:end) = [];
trl_n = trl_n(sel_idx);

chck = strcmp(ID2,ID(sel_idx)');
if any(chck==0)
    error('error number of files does not match');
end;

for jt = 1:length(sel_idx)
    id = [raw{sel_idx(jt),3},'.',num2str(raw{sel_idx(jt),4}),'.mat'];
    id2 = [ID2{jt},'.',num2str(trl_n(jt)),'.mat'];
    
    if ~strcmp(id,id2)
        error('file assignment does not match');
    end;
end;
% %%
% load('~/froux/chronset/thresholds/greedy_optim_NP_data_BCN_17-Sep-2015.mat');
% %tresh_params = cell(size(optim_data.hist_t,3),1);
% itresh = cell(size(optim_data.hist_t,3),1);
% for it = 1:size(optim_data.hist_t,3)
%     %tresh_params(it) = {squeeze(mean(mean(optim_data.hist_t(:,:,it),2),1))};
%     itresh(it) = {squeeze(mean(mean(optim_data.hist_t(:,:,it),2),1))};
% end;
%%
dat = cell(length(ID2),1);
parfor jt = 1:length(ID2)
    
    dum =load([path2featfiles,ID2{jt},'.',num2str(trl_n(jt)),'.mat']);
    dat{jt} = dum.savedata;
end;
clear dum;
%%
Y = mean(pred(sel_idx,:),2);

outL = find(sign(Y-2000)==1);
Y(outL) = [];
dat(outL) = [];

if length(Y) ~= length(dat)
    error('vactor must be same length');
end;
%%
itresh = cell(6,1);
itresh{1} = .1;%amplitude
itresh{2} = .9;%wiener entropy
itresh{3} = .1;%spectral change
itresh{4} = .1;% amplitude modulation
itresh{5} = .9;%frequency modulation
itresh{6} = .1;%goodness of pitch

isyl_t.t1 = 0.25;
isyl_t.t2 = 0.035;

inu = 4;

% number of permutations
Nper = 1000;

% histroy
hist_e = zeros(Ntrl,Nper,1);
hist_o = zeros(Ntrl,Nper,6);
hist_t = zeros(Ntrl,Nper,6);
test_e = zeros(Ntrl,Nper,1);
init_mle = zeros(Ntrl,Nper,1);

ei_training = zeros(Ntrl,1);
ei_test = zeros(Ntrl,1);

% freeze point step size
fp = 50;

%omega parameters
omega_orig1 = 0.05*ones(6,1);
omega_scale = 1;

%rate parameters
rdec = 0.5;
rinc = 1.1;

% training/test partitioning %
test_pct = 20;
%%
if f ==1
    parfor ot = 1:Ntrl
        
        %reset omega
        [omega] = omega_orig1;
        
        % median manual scores
        [mY1] = Y';
        
        %parameter space
        [param_space] =  repmat({[0 1]},[6 1]);
        
        % generate the partition of the training and test data
        idx = randperm(length(dat));
        idx = idx(1:round(length(dat)/100*test_pct));
        [test_idx] =  idx;
        [training_idx] = setdiff(1:length(dat),idx);
        
        %onsets for training data (initial guess treshold values)
        on1 = zeros(length(dat),1);
        for it = 1:length(dat)
            [on1(it)] =  detect_speech_on_and_offset_orig2(dat{it},[itresh' {0.035} {4}]');
        end;
        
        [mle1] = mle(on1(training_idx)-mY1(training_idx)');%mle for training data
        ei_training(ot) = mle1(2);
        
        [mle2] = mle(on1(test_idx)-mY1(test_idx)');% mle for test data
        ei_test(ot) = mle2(2);
        
        [test_ref] = mle2(2);
        
        [ci1] = mle1(2);% save the std of the training mle
        
        if ~isnan(ci1) && ~any(isnan(on1(training_idx)-mY1(training_idx)'))
            
            % reference parameters
            [ref_ci] = ci1;% mle for intial values
            [ref_tresh1] = itresh;% initial treshold values
            
            [d_e] = zeros(Nper,1);% track history of mle over iterations
            
            % loop over optimization-iterations
            for kt = 1:Nper
                s_t = tic;
                
                % randomly pick a feature
                sel = randperm(6);
                sel = sel(1);
                
                % update the treshold values
                temp1 = ref_tresh1;
                temp2 = ref_tresh1;
                
                % apply positive omega -increase temperature
                if sign((temp1{sel}+omega(sel))-max(param_space{sel}))==-1
                    temp1{sel} = temp1{sel}+omega(sel);
                end;
                
                %compute onsets
                on_a = zeros(length(dat),1);
                for it = 1:length(dat)
                    [on_a(it)] =  detect_speech_on_and_offset_orig2(dat{it},temp1);
                end;
                
                % estimate mle for positive temperature change
                ep = mle(on_a(training_idx)-mY1(training_idx)');%mle for training data
                ci2_a = ep(2);% save the std
                
                % apply negative omega - decrease temperature
                if sign((temp2{sel}-omega(sel))-min(param_space{sel}))==1
                    temp2{sel} = temp2{sel}-omega(sel);
                end;
                
                %compute onsets
                on_b = zeros(length(dat),1);
                for it = 1:length(dat)
                    [on_b(it)] =  detect_speech_on_and_offset_orig2(dat{it},temp2);
                end;
                
                % estimate mle for negative temperature change
                ep = mle(on_b(training_idx)-mY1(training_idx)');%training
                ci2_b = ep(2);% save the std
                
                % chose change direction
                if sign(ci2_a - ci2_b)==-1 % if positive temperature reduces mle
                    gd = ci2_a;%positive change
                    temp = temp1;
                    on1 = on_a;
                else
                    gd = ci2_b;%negative change
                    temp = temp2;
                    on1= on_b;
                end;
                
                if ~any(isnan(on1(training_idx)-mY1(training_idx)')) && sign(ref_ci-gd) == 1
                    ref_ci = gd;% set new reference for mle
                    ref_tresh1 = temp;% set ne reference for treshold
                    hist_e(ot,kt) = ref_ci;% record mle history
                    omega(sel) = omega(sel)*rinc;% apply rate increment to omega
                    
                    %compute onsets test
                    test_ep = mle(on1(test_idx)-mY1(test_idx)');%mle for test data
                    test_e(ot,kt) = test_ep(2);%keep track of mle for test data
                    test_ref = test_e(ot,kt);
                else
                    hist_e(ot,kt) = ref_ci;
                    omega(sel) = omega(sel)*rdec;% apply rate decrement to omega
                    test_e(ot,kt) = test_ref;% record mle history
                end;
                
                d_e(kt) =  ref_ci;
                
                %freeze point
                if mod(kt,fp) == 0
                    
                    dx = diff(d_e((kt-fp)+1:kt));
                    if max(dx) == 0 % if no change has happened over step size
                        omega = omega./omega_scale;% scale omega
                    end;
                    
                end;
                
                % keep track of omega
                hist_o(ot,kt,:) = omega;
                hist_t(ot,kt,:) = [ref_tresh1{:}];
                
            end;
        end;
    end;
else
     for ot = 1:Ntrl
        
        %reset omega
        [omega] = omega_orig1;
        
        % median manual scores
        [mY1] = Y';
        
        %parameter space
        [param_space] =  repmat({[0 1]},[6 1]);
        
        % generate the partition of the training and test data
        idx = randperm(length(dat));
        idx = idx(1:round(length(dat)/100*test_pct));
        [test_idx] =  idx;
        [training_idx] = setdiff(1:length(dat),idx);
        
        %onsets for training data (initial guess treshold values)
        on1 = zeros(length(dat),1);
        parfor it = 1:length(dat)
            [on1(it)] =  detect_speech_on_and_offset_orig2(dat{it},itresh);
        end;
        
        [mle1] = mle(on1(training_idx)-mY1(training_idx)');%mle for training data
        ei_training(ot) = mle1(2);
        
        [mle2] = mle(on1(test_idx)-mY1(test_idx)');% mle for test data
        ei_test(ot) = mle2(2);
        
        [test_ref] = mle2(2);
        
        [ci1] = mle1(2);% save the std of the training mle
        
        if ~isnan(ci1) && ~any(isnan(on1(training_idx)-mY1(training_idx)'))
            
            % reference parameters
            [ref_ci] = ci1;% mle for intial values
            [ref_tresh1] = itresh;% initial treshold values
            
            [d_e] = zeros(Nper,1);% track history of mle over iterations
            
            % loop over optimization-iterations
            for kt = 1:Nper
                s_t = tic;
                
                % randomly pick a feature
                sel = randperm(6);
                sel = sel(1);
                
                % update the treshold values
                temp1 = ref_tresh1;
                temp2 = ref_tresh1;
                
                % apply positive omega -increase temperature
                if sign((temp1{sel}+omega(sel))-max(param_space{sel}))==-1
                    temp1{sel} = temp1{sel}+omega(sel);
                end;
                
                %compute onsets
                on_a = zeros(length(dat),1);
                parfor it = 1:length(dat)
                    [on_a(it)] =  detect_speech_on_and_offset_orig2(dat{it},temp1);
                end;
                
                % estimate mle for positive temperature change
                ep = mle(on_a(training_idx)-mY1(training_idx)');%mle for training data
                ci2_a = ep(2);% save the std
                
                % apply negative omega - decrease temperature
                if sign((temp2{sel}-omega(sel))-min(param_space{sel}))==1
                    temp2{sel} = temp2{sel}-omega(sel);
                end;
                
                %compute onsets
                on_b = zeros(length(dat),1);
                parfor it = 1:length(dat)
                    [on_b(it)] =  detect_speech_on_and_offset_orig2(dat{it},temp2);
                end;
                
                % estimate mle for negative temperature change
                ep = mle(on_b(training_idx)-mY1(training_idx)');%training
                ci2_b = ep(2);% save the std
                
                % chose change direction
                if sign(ci2_a - ci2_b)==-1 % if positive temperature reduces mle
                    gd = ci2_a;%positive change
                    temp = temp1;
                    on1 = on_a;
                else
                    gd = ci2_b;%negative change
                    temp = temp2;
                    on1= on_b;
                end;
                
                if ~any(isnan(on1(training_idx)-mY1(training_idx)')) && sign(ref_ci-gd) == 1
                    ref_ci = gd;% set new reference for mle
                    ref_tresh1 = temp;% set ne reference for treshold
                    hist_e(ot,kt) = ref_ci;% record mle history
                    omega(sel) = omega(sel)*rinc;% apply rate increment to omega
                    
                    %compute onsets test
                    test_ep = mle(on1(test_idx)-mY1(test_idx)');%mle for test data
                    test_e(ot,kt) = test_ep(2);%keep track of mle for test data
                    test_ref = test_e(ot,kt);
                else
                    hist_e(ot,kt) = ref_ci;
                    omega(sel) = omega(sel)*rdec;% apply rate decrement to omega
                    test_e(ot,kt) = test_ref;% record mle history
                end;
                
                d_e(kt) =  ref_ci;
                
                %freeze point
                if mod(kt,fp) == 0
                    
                    dx = diff(d_e((kt-fp)+1:kt));
                    if max(dx) == 0 % if no change has happened over step size
                        omega = omega./omega_scale;% scale omega
                    end;
                    
                end;
                
                % keep track of omega
                hist_o(ot,kt,:) = omega;
                hist_t(ot,kt,:) = [ref_tresh1{:}];
            end;
        end;
    end;
end;
%%
optim_data = struct;

optim_data.fp = fp;
optim_data.omega_orig = omega_orig1;
optim_data.omega_scale = omega_scale;
optim_data.rdec = rdec;
optim_data.rinc = rinc;
optim_data.test_pct =test_pct;
optim_data.Nruns = Ntrl;
optim_data.Niter = Nper;

% history of initial threshold values
optim_data.itresh = itresh;
% history of mle for training data
optim_data.hist_e = hist_e;clear hist_e;
% history of mle for test data
optim_data.test_e = test_e;clear test_e;
% history of threshold
optim_data.hist_t = hist_t;clear hist_t;
% history of omega
optim_data.hist_o = hist_o;clear hist_o;
% initial mle for training data
optim_data.ei_training = ei_training;clear ei_training;
% initial mle for test data
optim_data.ei_test = ei_test;clear ei_test;

path2files = '/bcbl/home/home_a-f/froux/chronset/thresholds/';
ct = clock;
ct = ct(4:6);
ct = num2str(round(ct));
save([path2files,'greedy_optim_SayWhen_',date,'_',c(1),'-',c(2),'-x',c(3),'.mat'],'optim_data');
%%
matlabpool('close');











