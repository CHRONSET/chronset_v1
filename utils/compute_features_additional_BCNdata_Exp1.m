%%
if matlabpool('size')==0
    matlabpool 64;%ips_base;
end;
%%
restoredefaultpath;
addpath(genpath('/bcbl/home/home_a-f/froux/chronset/'));
%%
[path2files] = '~froux/froux/chronset/data/BCN/additional_data/raw_data/';
[files] = dir([path2files,'*.WAV']);
[save_path] = '~froux/froux/chronset/data/BCN/additional_data/raw_data/';
%%
ts = tic;
parfor jt = 1:length(files)

    feat_data = struct;
    [feat_data.wav,feat_data.FS] = wavread2([path2files,files(jt).name]);
    feat_data.wav = feat_data.wav(:,1);
    
    %[feat_data2] = compute_feat_data(feat_data);
    [feat_data2] = compute_feat_data(feat_data);
    
    dat = struct;
    dat.features = feat_data2.features;
    dat.t = feat_data2.t;
    dat.syl_seg = feat_data2.syl_seg;
    
    savename = [files(jt).name(1:end-3),'mat'];
    savename(regexp(savename,' ')) = [];
    
    save_parfor(save_path,savename,dat);
    
end;
toc(ts);
clear feat_data2 dat;
%%
matlabpool('close');
%%
exit;