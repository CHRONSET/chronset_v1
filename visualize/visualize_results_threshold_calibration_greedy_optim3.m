%%
restoredefaultpath;
addpath(genpath('/bcbl/home/home_a-f/froux/chronset/'));
%%
p2df = '/bcbl/home/home_a-f/froux/chronset/thresholds/';
fn ='greedy_optim_SayWhen_18-Sep-2015.mat';
load([p2df,fn]);
%%
optim_data.test_e(optim_data.test_e==0) = NaN;
ind = find(optim_data.test_e == min(min(optim_data.test_e)));
[i1,i2] = ind2sub(size(optim_data.test_e),ind);
%%
i1 = min(unique(i1));
i2 = min(unique(i2));

[i1 i2]
%%
x = squeeze([optim_data.hist_t(i1,i2,:)]);

t = cell(length(x),1);
for it = 1:length(x)
    t{it} = x(it);
end;
%%
figure;
subplot(2,2,1);
hold on;
plot([optim_data.itresh{:}],'rs-');
plot(x,'bo--');
legend('training','test');
title('threshold values');
subplot(2,2,2);
plot(squeeze(optim_data.hist_o(i1,:,:)));
title('\Omega');

subplot(223);
hold on;
plot(1:size(optim_data.hist_e,2),optim_data.hist_e,'b-');
plot(1:size(optim_data.hist_e,2),optim_data.hist_e(i1,:),'k-','LineWidth',3);
title('training error');

subplot(224);
hold on;
plot(1:size(optim_data.test_e,2),optim_data.test_e,'r-');
plot(1:size(optim_data.test_e,2),optim_data.test_e(i1,:),'k-','LineWidth',3);
title('threshold values');
title('test error');

%%
figure;
for it = 1:size(optim_data.hist_t,3)
    subplot(5,2,it);
    hold on;
    plot(1:size(optim_data.hist_t,2),squeeze([optim_data.hist_t(:,:,it)]),'Color',[.75 .75 .75]);
    plot(1:size(optim_data.hist_t,2),squeeze([optim_data.hist_t(i1,:,it)]),'Color',[0 0 0],'LineWidth',3);

end;
% %%
% readout_manual_ratings_additionalBCNdata;
% %%
% if matlabpool('size')==0
%     matlabpool local;
% end;
% 
% p2df = '/bcbl/home/home_a-f/froux/chronset/data/BCN/additional_data/feature_data/Jun2016/';
% 
% idx1 = optim_data.hist_t(i1,i2,9);
% idx2 = optim_data.hist_t(i1,i2,10);
% 
% dat = cell(length(fID),1);
% k = 0;
% parfor it = 1:length(fID)
%         
%     chck = dir([p2df,[fID{it}(1:end-3),'_*Jun-2016.mat']]);
%     dum = load([p2df,chck.name]);
%     
%     dat{it} = dum.savedata{idx1,idx2};
%     dat{it}.id = fID{it};
%     
% end;
% clear dum;
% %%
% on = zeros(length(dat),1);
% for it = 1:length(dat)
%     [on(it)] =  detect_speech_on_and_offset_orig2(dat{it},t);
% end;
% %%
% figure;
% plot(on,man_RT,'bo');
% %%
% fid = fopen('data_files_list.txt','w+');
% for it = 1:length(fID)
%     fprintf(fid,'%s',fID{it});
%     fprintf(fid,'\n');
% end;
