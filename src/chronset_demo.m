function chronset_demo(input_file_name,varargin)
if nargin == 0
    input_file_name = '/bcbl/home/home_a-f/froux/chronset/data/demo/1P8_91235211.WAV';
end;
%%
restoredefaultpath;
addpath(genpath('/bcbl/home/home_a-f/froux/chronset/'));
%%
if isempty(varargin)
    [path2file,name,~] = fileparts(input_file_name);
    savepath = [path2file,'/'];
end;
%%
in = struct;
[in.wav,in.FS] = wavread2(input_file_name);
%%
load('/bcbl/home/home_a-f/froux/chronset/thresholds/greedy_optim_NP_data_BCN_03-Sep-2015.mat');
[i1,i2] = find(optim_data.hist_e == min(min(optim_data.hist_e)));
i1 = min(unique(i1));
i2 = min(unique(i2));
%%
tresh{1} = squeeze(optim_data.hist_t(i1,i2,1));
tresh{2} = squeeze(optim_data.hist_t(i1,i2,2));
tresh{3} = squeeze(optim_data.hist_t(i1,i2,3));
tresh{4} = squeeze(optim_data.hist_t(i1,i2,4));
tresh{5} = squeeze(optim_data.hist_t(i1,i2,5));
tresh{6} = squeeze(optim_data.hist_t(i1,i2,6));
%%
[feat_data] = compute_feat_data(in);
% %%
% amp = feat_data.features{1};
% we = feat_data.features{2};
% dfdt = (feat_data.features{3});
% am = feat_data.features{4};
% fm = feat_data.features{5};
% goP = feat_data.features{6};
% 
% X = [amp we dfdt am fm goP];
% 
% [u,s,v] = svd(X);
% 
% dum = u(:,1)*s(1,1)*v(:,1)';
%%
[on] = detect_speech_on_and_offset(feat_data,tresh);
%%
fid = fopen([savepath,name,'_onset.txt'],'w+');
fprintf(fid,'%s',name);
fprintf(fid,'\t');
fprintf(fid,'%s',num2str(on));
fprintf(fid,'\n\r');
