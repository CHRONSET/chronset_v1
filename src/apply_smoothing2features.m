function [sm_feat_data]= apply_smoothing2features(feat_data,kernel_length,kernel_type)

dum = [feat_data(1)*ones(1,1000) feat_data' feat_data(end)*ones(1,1000)];

dum = eval(['conv(dum,',kernel_type,'(',num2str(kernel_length),'))./sum(',kernel_type,'(',num2str(kernel_length),'))']);

dum = dum';

[sm_feat_data]= dum((1000+(kernel_length/2)+1):end-(1000+(kernel_length/2)-1));

%TO DO
% the kernel needs to be normalized to 1