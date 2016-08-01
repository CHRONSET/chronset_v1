%% 
function [thresh] = chronset_extract_thresholds(optim_data)

[i1,i2] = find(optim_data.hist_e == min(min(optim_data.hist_e)));

i1 = unique(i1);
i2 = unique(i2);

i1 = i1(1);
i2 = i2(1);

fprintf(['the selected thresholds have a training error of: ',num2str(optim_data.hist_e(i1,i2)),'ms\n']);
fprintf(['the selected thresholds have a testing error of: ',num2str(optim_data.test_e(i1,i2)),'ms\n']);

fprintf(['loading feature thresholds\n']);
thresh = cell(size(optim_data.hist_t,3),1);
for it = 1:length(thresh)
    thresh{it} = squeeze(optim_data.hist_t(i1,i2,it));
end;
fprintf(['\n']);
fprintf(['\n']);