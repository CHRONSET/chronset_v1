addpath(genpath('~/froux/chronset/src/'));
%%
load('/bcbl/home/home_a-f/froux/chronset/thresholds/greedy_optim_SayWhen_18-Sep-2015');

[idx1,idx2] = find(optim_data.hist_e == min(min(optim_data.hist_e)));

[tparams] = squeeze(optim_data.hist_t(median(idx1),floor(median(idx2)),:));

figure;
subplot(131);
plot(squeeze(optim_data.hist_o(1,:,:)));
subplot(132);
hold on;
plot(mean(optim_data.hist_e(1,:),1));
plot(mean(optim_data.test_e(1,:),1),'r');
subplot(133);
hold on;
%plot(trck_ml1);
%plot(trck_ml2,'r');
%%
p = [[optim_data.itresh{:}]' tparams];
p = [p diff(p,[],2)]
%%
trsh = cell(length(tparams),1);
for it = 1:length(tparams)
    trsh{it}= tparams(it);
end;

on = zeros(length(dat),1);
parfor it = 1:length(dat)

    %[on(it)] = detect_speech_on_and_offset(dat{it},trsh,[]);
    [on(it)] = detect_speech_on_and_offset_orig2(dat{it},trsh);

end;
%%
mY1 = mean(Y,2);

figure;plot(on,mY1,'b.');
[b,~,~,~,rsq] = regress(mY1,[ones(size(on)) on]);
rsq(1)