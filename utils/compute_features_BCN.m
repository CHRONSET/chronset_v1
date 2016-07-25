%%
if matlabpool('size')==0
    matlabpool 96;
end;
%%
restoredefaultpath;
addpath(genpath('/bcbl/home/home_a-f/froux/chronset/'));
%%
[path2files] = '/bcbl/home/home_a-f/froux/chronset/data/BCN/raw_data/';
[files] = dir([path2files,'*.WAV']);
[save_path] = '/bcbl/home/home_a-f/froux/chronset/data/BCN/feature_data/';
% %%
% chck = dir([save_path,'sel_idx_*.mat']);
% 
% if isempty(chck)
%     sel_idx = randperm(length(files));
%     sel_idx = sel_idx(1:400);
% 
%     save([save_path,'sel_idx_',date,'.mat'],'sel_idx');
% else
%     load([save_path,chck.name]);
% end;
% 
% files = files(sel_idx);
%%
p1 = [.004:.003:.016];
p2 = [250:250:1250];
%%
ts = tic;
%for jt = 1;
parfor jt = 1:length(files)
    
    feat_data = struct;
    [feat_data.wav,feat_data.FS] = wavread2([path2files,files(jt).name]);
    feat_data.wav = feat_data.wav(:,1);
    
    optp =[];
    dat = cell(length(p1),length(p2));
    for nt = 1:length(p1)
        
        optp.p1 = p1(nt);
        
        for ot = 1:length(p2)
            
            optp.p2 = p2(ot);
            
            %[feat_data2] = compute_feat_data(feat_data);
            [feat_data2] = compute_feat_data(optp,feat_data);
            

            dat{nt,ot}.features = feat_data2.features;
            dat{nt,ot}.t = feat_data2.t;
            dat{nt,ot}.syl_seg = feat_data2.syl_seg;
            
        end;
    end;
   
    savename = [files(jt).name(1:end-4),'_',date,'.mat'];
    savename(regexp(savename,' ')) = [];
    
    save_parfor(save_path,savename,dat);

end;
toc(ts);
clear feat_data2 dat;
%%
matlabpool('close');
exit;
