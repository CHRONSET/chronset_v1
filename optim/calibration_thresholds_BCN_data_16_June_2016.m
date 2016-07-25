function calibration_thresholds_BCN_data_16_June_2016(f,Nper,nworker,lf)
%%
if nargin == 0
    Nper = 500;% number of optim iterations
    f = 1;
    nworker = 32;
    lf = 8;
end;
%%
restoredefaultpath;
addpath(genpath('/bcbl/home/home_a-f/froux/chronset/'));
% % %% init workspace
% clear;
% clc;
% close all;
%%
if f == 1
    if matlabpool('size') ==0
        matlabpool(nworker);% number of training vs test partitions
    end;
else
    if matlabpool('size') ==0
        matlabpool local;
    end;
end;
%% read in the manual scores
path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/';

[data,txt,raw] = xlsread([path2files,'manual_ratings_BCN.xls']);

% get rid of header
txt(1,:) = [];
raw(1,:) = [];

% manual onsets
Y = cell2mat(raw(:,9:10));


[pc] = txt(:,1); % participant ID
[fc] = raw(:,4); % file ID

% convert cell to double
fc2 = zeros(length(fc)-1,1);
parfor it = 1:length(fc)
    fc2(it) = fc{it};
end;

if length(pc) ~= length(fc2) || length(pc) ~= length(Y) || length(fc) ~=length(Y)
    error('wrong number of files');
end;
%% generate ID information
path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/raw_data/';
ID = cell(length(fc2),1);
parfor it = 1:length(fc2)
    ID{it} = [pc{it},'_',num2str(fc2(it)),'.WAV'];
end;


[sel_idx] = zeros(length(ID),1);
parfor it = 1:length(ID)
    file_name = dir([path2files,ID{it}]);
    if ~isempty(file_name)
        sel_idx(it) = it;
    end;
end;
[sel_idx] = sel_idx(sel_idx ~= 0);

if (length(sel_idx) ~= length(pc)) || (length(sel_idx) ~= length(fc2)) || length(sel_idx) ~= length(Y)
    error('Wrong number of files');
end;
%% generate ID information
path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/feature_data/Jun2016/';
ID2 = cell(length(fc2),1);
parfor it = 1:length(fc2)
     fn = dir([path2files,pc{it},'_',num2str(fc2(it)),'*-Jun-2016.mat']);
     ID2{it} = fn.name;
end;


[sel_idx] = zeros(length(ID2),1);
parfor it = 1:length(ID2)
    file_name = dir([path2files,ID2{it}]);
    if ~isempty(file_name)
        sel_idx(it) = it;
    end;
end;
[sel_idx] = sel_idx(sel_idx ~= 0);

if (length(sel_idx) ~= length(pc)) || (length(sel_idx) ~= length(fc2)) || (length(sel_idx) ~= length(Y))
    error('Wrong number of files');
end;
%%
chck = zeros(length(ID2),1);
parfor it = 1:length(ID2)
    chck(it) = strcmp(ID{it}(1:end-4),ID2{it}(1:min(regexp(ID2{it}(1:end-4),'-'))-4));
end;
if sum(chck) ~= length(ID2)
    error('number of files does not match');
end;
%%
path2files = '/bcbl/home/home_a-f/froux/chronset/data/BCN/feature_data/Jun2016/';

dat = cell(length(ID2),2,2);
k = 0;
parfor it = 1:length(ID2)
    
    dum = load([path2files,ID2{it}]);
    
    dat(it,:,:) = [dum.savedata(3,3) dum.savedata(3,4);dum.savedata(4,3) dum.savedata(4,4)];
    
end;
clear dum;
%%
% number of partitions training vs test
if f ==1
    Ntrl = matlabpool('size')*lf;
    
    %starting values for thresholds
    itresh = cell(10,1);
    itresh{1} = .1;%amplitude
    itresh{2} = .9;%wiener entropy
    itresh{3} = .1;%spectral change
    itresh{4} = .1;% amplitude modulation
    itresh{5} = .9;%frequency modulation
    itresh{6} = .1;%goodness of pitch
    
%     %added June 13 2016
    itresh{7}  = .035;%syl length
    itresh{8}  = 4;%number of simult active features
    
    itresh{9} = 1;%sliding window size
    itresh{10} = 1;%smoothing
    
    nfeat = length(itresh);
    
    % histroy
    hist_e = zeros(Ntrl,Nper,1);
    hist_o = zeros(Ntrl,Nper,nfeat);
    hist_t = zeros(Ntrl,Nper,nfeat);
    test_e = zeros(Ntrl,Nper,1);
    init_mle = zeros(Ntrl,Nper,1);
    
    trck_ml1 = zeros(Ntrl,Nper,1);
    trck_ml2 = zeros(Ntrl,Nper,1);
    
    ei_training = zeros(Ntrl,1);
    ei_test = zeros(Ntrl,1);
    
    % freeze point step size
    fp = 50;
    
    %omega parameters
    omega_orig1 = 0.05*ones(6,1);
    if sign(length(itresh)-6) ==1
        omega_orig1(7:10) = [itresh{7:10}];
    end;
    
    omega_scale = 1;
    
    %rate parameters
    rdec = 0.5;
    rinc = 1.1;
    
    % training/test partitioning %
    test_pct = 20;
    
    Y = Y(sel_idx,:);
    %%
    s_t = tic;
    parfor ot = 1:Ntrl

        %reset omega
            [omega] = omega_orig1;
            
            % median manual scores
            [mY1] = median(Y,2)';
            
            %parameter space
            [param_space]   =  repmat({[0 1]},[6 1]);
            
            %added June 13 2016
            if sign(length(itresh)-6) ==1
                
                param_space{7}  = [.005 .05];
                param_space{8}  = [3    5];
                
                param_space{9}  = [1   2];
                param_space{10}  = [1   2];
            end;
            nchck = 0;
            % generate the partition of the training and test data
            idx = randperm(length(dat));
            idx = idx(1:round(length(dat)/100*test_pct));
            [test_idx] =  idx;
            [training_idx] = setdiff(1:length(dat),idx);
            
            %onsets for training data (initial guess treshold values)
            on1 = zeros(length(dat),1);
            if sign(nfeat-6)==1
                for it = 1:length(dat)
                    [on1(it)] =  detect_speech_on_and_offset_orig2(dat{it,itresh{9},itresh{10}},itresh);
                end;
            else
                for it = 1:length(dat)
                    [on1(it)] =  detect_speech_on_and_offset_orig2(dat{it,1,1},[itresh{:} {.035} {4}]);
                end;
            end;
            
            %         if any(isnan(on1))
            %             error('initial seed must be integer');
            %         end;
            
            
            %[mle1] = mle(on1(training_idx)-mY1(training_idx)');%mle for training data
            C = mY1(training_idx)';
            X = on1(training_idx);
            [~,~,r1,~,~] = regress(C,[ones(size(X)) X]);
            [mle1] = mle(r1);
            ei_training(ot) = mle1(2);
            
            %[mle2] = mle(on1(test_idx)-mY1(test_idx)');% mle for test data
            C = mY1(test_idx)';
            X = on1(test_idx);
            [~,~,r2,~,~] = regress(C,[ones(size(X)) X]);
            [mle2] = mle(r2);
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
                    
                    % randomly pick a feature                             
                    sel = randperm(nfeat);
                    sel = sel(1);                    
                    
                    % update the treshold values
                    temp1 = ref_tresh1;
                    temp2 = ref_tresh1;
                    
                    % apply positive omega -increase temperature
                    if ismember(sel,1:6)
                        if sign((temp1{sel}+omega(sel))-max(param_space{sel}))==-1
                            temp1{sel} = temp1{sel}+omega(sel);
                        end;
                    end;
                    
                    if sign(length(itresh)-6) ==1
                        
                        if (sel == 7)
                            if sign((temp1{sel}+.001)-max(param_space{sel}))==-1
                                temp1{sel} = temp1{sel}+.001;
                            end;
                        end;
                        
                        if ismember(sel,8:10)
                            if sign((temp1{sel}+1)-max(param_space{sel}))==-1 || sign((temp1{sel}+1)-max(param_space{sel}))== 0
                                temp1{sel} = temp1{sel}+1;
                            end;
                        end;
                    end;
                    
                    %compute onsets
                    on_a = zeros(length(dat),1);
                    if sign(nfeat-6)==1
                        for it = 1:length(dat)
                            [on_a(it)] =  detect_speech_on_and_offset_orig2(dat{it,temp1{9},temp1{10}},temp1);
                        end;
                    else
                        for it = 1:length(dat)
                            [on_a(it)] =  detect_speech_on_and_offset_orig2(dat{it,1,1},[temp1{:} {.035} {4}]);
                        end;
                    end;
                    
                    %                 if any(isnan(on_a))
                    %                     error('onset a must be integer');
                    %                 end;
                    
                    % estimate mle for positive temperature change
                    %ep = mle(on_a(training_idx)-mY1(training_idx)');%mle for training data
                    [~,~,r,~,~] = regress(mY1(training_idx)',[ones(size(on_a(training_idx))) on_a(training_idx)]);
                    ep = mle(r);
                    ci2_a = ep(2);% save the std
                    
                    trck_ml1(ot,kt) = ep(2);
                    
                    % modified June 13 2016
                    % apply negative omega - decrease temperature
                    if ismember(sel,1:6)
                        if sign((temp2{sel}-omega(sel))-min(param_space{sel}))==1
                            temp2{sel} = temp2{sel}-omega(sel);
                        end;
                    end;
                    
                    if sign(length(itresh)-6) ==1     
                        
                        if (sel == 7)
                            if sign((temp2{sel}-.001)-min(param_space{sel}))==1
                                temp2{sel} = temp2{sel}-.001;
                            end;
                        end;
                        
                        if ismember(sel,8:10)
                            if (sign((temp2{sel}-1)-min(param_space{sel}))==1) || (sign((temp2{sel}-1)-min(param_space{sel}))== 0)
                                temp2{sel} = temp2{sel}-1;
                            end;
                        end;
                        
                    end;
                    
                    %compute onsets
                    on_b = zeros(length(dat),1);
                    if sign(nfeat-6)==1
                        for it = 1:length(dat)
                            [on_b(it)] =  detect_speech_on_and_offset_orig2(dat{it,temp2{9},temp2{10}},temp2);
                        end;
                    else
                        for it = 1:length(dat)
                            [on_b(it)] =  detect_speech_on_and_offset_orig2(dat{it,1,1},[temp2{:} {.035} {4}]);
                        end;
                    end;
                    
                    %                 if any(isnan(on_b))
                    %                     error('onset b seed must be integer');
                    %                 end;
                    
                    % estimate mle for negative temperature change
                    %ep = mle(on_b(training_idx)-mY1(training_idx)');%training
                    [~,~,r,~,~] = regress(mY1(training_idx)',[ones(size(on_b(training_idx))) on_b(training_idx)]);
                    ep = mle(r);
                    ci2_b = ep(2);% save the std
                    
                    trck_ml2(ot,kt) = ep(2);
                    
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
                    
                    x = [0 1];draw = randperm(2);x = x(draw(1));%draw random num
                    
                    if ( ~any(isnan(on1(training_idx)-mY1(training_idx)')) && sign(ref_ci-gd) == 1 ) || ( ( ~any(isnan(on1(training_idx)-mY1(training_idx)')) && sign(ref_ci-gd) == 0 ) && ( x==1 ) )
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
    Ntrl = 1;
    
    %starting values for thresholds
    itresh = cell(6,1);
    itresh{1} = .1;%amplitude
    itresh{2} = .9;%wiener entropy
    itresh{3} = .1;%spectral change
    itresh{4} = .1;% amplitude modulation
    itresh{5} = .9;%frequency modulation
    itresh{6} = .1;%goodness of pitch
    
    % histroy
    hist_e = zeros(Ntrl,Nper,1);
    hist_o = zeros(Ntrl,Nper,6);
    hist_t = zeros(Ntrl,Nper,6);
    test_e = zeros(Ntrl,Nper,1);
    init_mle = zeros(Ntrl,Nper,1);
    
    trck_ml1 = zeros(Ntrl,Nper,1);
    trck_ml2 = zeros(Ntrl,Nper,1);
    
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
    
    %Y = Y(sel_idx,:);
    %%
    s_t = tic;
    for ot = 1:Ntrl
        %reset omega
        [omega] = omega_orig1;
        
        % median manual scores
        [mY1] = mean(Y,2)';
        
        %parameter space
        [param_space] =  repmat({[0 1]},[6 1]);
        
        nchck = 0;
        % generate the partition of the training and test data
        idx = randperm(length(dat));
        idx = idx(1:round(length(dat)/100*test_pct));
        [test_idx] =  idx;
        [training_idx] = setdiff(1:length(dat),idx);
        
        %onsets for training data (initial guess treshold values)
        on1 = zeros(length(dat),1);
        parfor it = 1:length(dat)
            [on1(it)] =  detect_speech_on_and_offset(dat{it},itresh,mY1(it));
            %[on1(it)] =  detect_speech_on_and_offset_orig2(dat{it},itresh);
        end;
        
        %         if any(isnan(on1))
        %             error('initial seed must be integer');
        %         end;
        
        %[mle1] = mle(on1(training_idx)-mY1(training_idx)');%mle for training data
        C = mY1(training_idx)';
        X = on1(training_idx);
        [~,~,r1,~,~] = regress(C,[ones(size(X)) X]);
        [mle1] = mle(r1);
        ei_training(ot) = mle1(2);
        
        %[mle2] = mle(on1(test_idx)-mY1(test_idx)');% mle for test data
        C = mY1(test_idx)';
        X = on1(test_idx);
        [~,~,r2,~,~] = regress(C,[ones(size(X)) X]);
        [mle2] = mle(r2);
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
                    [on_a(it)] =  detect_speech_on_and_offset(dat{it},temp1,mY1(it));
                    %[on_a(it)] =  detect_speech_on_and_offset_orig2(dat{it},temp1);
                end;
                
                %                 if any(isnan(on_a))
                %                     error('onset a must be integer');
                %                 end;
                
                % estimate mle for positive temperature change
                %ep = mle(on_a(training_idx)-mY1(training_idx)');%mle for training data
                [~,~,r,~,~] = regress(mY1(training_idx)',[ones(size(on_a(training_idx))) on_a(training_idx)]);
                ep = mle(r);
                ci2_a = ep(2);% save the std
                
                trck_ml1(ot,kt) = ep(2);
                
                % apply negative omega - decrease temperature
                if sign((temp2{sel}-omega(sel))-min(param_space{sel}))==1
                    temp2{sel} = temp2{sel}-omega(sel);
                end;
                
                %compute onsets
                on_b = zeros(length(dat),1);
                parfor it = 1:length(dat)
                    [on_b(it)] =  detect_speech_on_and_offset(dat{it},temp2,mY1(it));
                    %[on_b(it)] =  detect_speech_on_and_offset_orig2(dat{it},temp2);
                    
                end;
                
                %                 if any(isnan(on_b))
                %                     error('onset b seed must be integer');
                %                 end;
                
                % estimate mle for negative temperature change
                %ep = mle(on_b(training_idx)-mY1(training_idx)');%training
                [~,~,r,~,~] = regress(mY1(training_idx)',[ones(size(on_b(training_idx))) on_b(training_idx)]);
                ep = mle(r);
                ci2_b = ep(2);% save the std
                
                trck_ml2(ot,kt) = ep(2);
                
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
%%
path2files = '/bcbl/home/home_a-f/froux/chronset/thresholds/';
save([path2files,'greedy_optim_NP_data_BCN_',date,'.mat'],'optim_data');
%%
matlabpool('close');








