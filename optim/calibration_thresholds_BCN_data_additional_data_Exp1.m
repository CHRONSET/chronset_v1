%%
restoredefaultpath;
addpath(genpath('/bcbl/home/home_a-f/froux/chronset/'));
%% init workspace
clear;
clc;
close all;
f = 1;
if f == 1
    if matlabpool('size') ==0
        matlabpool 100;
    end;
else    
    if matlabpool('size') ==0
        matlabpool local;
    end;
end;
%%
readout_manual_ratings_additionalBCNdata;
%%
p2df = '/bcbl/home/home_a-f/froux/chronset/data/BCN/additional_data/feature_data/';

dat = cell(length(fID),1);
k = 0;
parfor it = 1:length(fID)
        
    dum = load([p2df,[fID{it}(1:end-3),'mat']]);
    dat{it} = dum.savedata;
    dat{it}.id = fID{it};
    
end;
clear dum;
%%
% number of partitions training vs test
if f ==1
    Ntrl = matlabpool('size')*1;
    % number of optim iterations
    Nper = 1000;
    
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
    
    Y = man_RT;
    %%
    s_t = tic;
    parfor ot = 1:Ntrl
        %reset omega
        [omega] = omega_orig1;       
        
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
        for it = 1:length(dat)
            [on1(it)] =  detect_speech_on_and_offset_orig2(dat{it},itresh);
        end;
        
%         if any(isnan(on1))
%             error('initial seed must be integer');
%         end;
        

        %[mle1] = mle(on1(training_idx)-Y(training_idx));%mle for training data
        C = Y(training_idx);
        X = on1(training_idx);
        [~,~,r1,~,~] = regress(C,[ones(size(X)) X]);
        [mle1] = mle(r1);
        ei_training(ot) = mle1(2);
        
        %[mle2] = mle(on1(test_idx)-Y(test_idx));% mle for test data
        C = Y(test_idx);
        X = on1(test_idx);
        [~,~,r2,~,~] = regress(C,[ones(size(X)) X]);
        [mle2] = mle(r2);
        ei_test(ot) = mle2(2);
        
        [test_ref] = mle2(2);
        
        [ci1] = mle1(2);% save the std of the training mle
        
        if ~isnan(ci1) && ~any(isnan(on1(training_idx)-Y(training_idx)))
            
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
                for it = 1:length(dat)
                    [on_a(it)] =  detect_speech_on_and_offset_orig2(dat{it},temp1);
                end;
                
%                 if any(isnan(on_a))
%                     error('onset a must be integer');
%                 end;
                
                % estimate mle for positive temperature change
                %ep = mle(on_a(training_idx)-Y(training_idx));%mle for training data
                [~,~,r,~,~] = regress(Y(training_idx),[ones(size(on_a(training_idx))) on_a(training_idx)]);
                ep = mle(r);
                ci2_a = ep(2);% save the std
                
                trck_ml1(ot,kt) = ep(2);
                
                % apply negative omega - decrease temperature
                if sign((temp2{sel}-omega(sel))-min(param_space{sel}))==1
                    temp2{sel} = temp2{sel}-omega(sel);
                end;
                
                %compute onsets
                on_b = zeros(length(dat),1);
                for it = 1:length(dat)
                    [on_b(it)] =  detect_speech_on_and_offset_orig2(dat{it},temp2);
                end;
                
%                 if any(isnan(on_b))
%                     error('onset b seed must be integer');
%                 end;
                
                % estimate mle for negative temperature change
                %ep = mle(on_b(training_idx)-Y(training_idx));%training
                [~,~,r,~,~] = regress(Y(training_idx),[ones(size(on_b(training_idx))) on_b(training_idx)]);
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
                
                if ~any(isnan(on1(training_idx)-Y(training_idx))) && sign(ref_ci-gd) == 1
                    ref_ci = gd;% set new reference for mle
                    ref_tresh1 = temp;% set ne reference for treshold
                    hist_e(ot,kt) = ref_ci;% record mle history
                    omega(sel) = omega(sel)*rinc;% apply rate increment to omega
                    
                    %compute onsets test
                    test_ep = mle(on1(test_idx)-Y(test_idx));%mle for test data
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
    % number of optim iterations
    Nper = 500;
    
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
    
    Y = man_RT;
    %%
    s_t = tic;
    for ot = 1:Ntrl
        %reset omega
        [omega] = omega_orig1;       
        
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
            [on1(it)] =  detect_speech_on_and_offset(dat{it},itresh,Y(it));
            %[on1(it)] =  detect_speech_on_and_offset_orig2(dat{it},itresh);
        end;
        
%         if any(isnan(on1))
%             error('initial seed must be integer');
%         end;
        
        %[mle1] = mle(on1(training_idx)-Y(training_idx));%mle for training data
        C = Y(training_idx);
        X = on1(training_idx);
        [~,~,r1,~,~] = regress(C,[ones(size(X)) X]);
        [mle1] = mle(r1);
        ei_training(ot) = mle1(2);
        
        %[mle2] = mle(on1(test_idx)-Y(test_idx));% mle for test data
        C = Y(test_idx);
        X = on1(test_idx);
        [~,~,r2,~,~] = regress(C,[ones(size(X)) X]);
        [mle2] = mle(r2);
        ei_test(ot) = mle2(2);
                
        [test_ref] = mle2(2);
        
        [ci1] = mle1(2);% save the std of the training mle
        
        if ~isnan(ci1) && ~any(isnan(on1(training_idx)-Y(training_idx)))
            
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
                    [on_a(it)] =  detect_speech_on_and_offset(dat{it},temp1,Y(it));
                    %[on_a(it)] =  detect_speech_on_and_offset_orig2(dat{it},temp1);
                end;
                
%                 if any(isnan(on_a))
%                     error('onset a must be integer');
%                 end;
                
                % estimate mle for positive temperature change
                %ep = mle(on_a(training_idx)-Y(training_idx));%mle for training data
                [~,~,r,~,~] = regress(Y(training_idx),[ones(size(on_a(training_idx))) on_a(training_idx)]);
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
                    [on_b(it)] =  detect_speech_on_and_offset(dat{it},temp2,Y(it));
                    %[on_b(it)] =  detect_speech_on_and_offset_orig2(dat{it},temp2);
                    
                end;
                
%                 if any(isnan(on_b))
%                     error('onset b seed must be integer');
%                 end;
                
                % estimate mle for negative temperature change
                %ep = mle(on_b(training_idx)-Y(training_idx));%training
                [~,~,r,~,~] = regress(Y(training_idx),[ones(size(on_b(training_idx))) on_b(training_idx)]);
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
                
                if ~any(isnan(on1(training_idx)-Y(training_idx))) && sign(ref_ci-gd) == 1
                    ref_ci = gd;% set new reference for mle
                    ref_tresh1 = temp;% set ne reference for treshold
                    hist_e(ot,kt) = ref_ci;% record mle history
                    omega(sel) = omega(sel)*rinc;% apply rate increment to omega
                    
                    %compute onsets test
                    test_ep = mle(on1(test_idx)-Y(test_idx));%mle for test data
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
save([path2files,'greedy_optim_NP_data_BCN_additional_data_Exp1',date,'.mat'],'optim_data');
%%
matlabpool('close');