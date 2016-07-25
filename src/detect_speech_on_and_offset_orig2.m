function [on,off,feat_data] =  detect_speech_on_and_offset_orig2(feat_data,tresh)
% %%
% tresh = cell(length(feat_data.features),1);
% tresh{1} = 5e-3;%amplitude
% tresh{2} = .55;%wiener entropy
% tresh{3} = 5e-3;%spectral change
% tresh{4} = 5e-3;% amplitude modulation
% tresh{5} = .6;%frequency modulation
% tresh{6} = .25;%goodness of pitch

% new features for optim add June 13 2016
% tresh{7}  word length
% thresh{8} sylable length
% thres{9} number of simultaneous features
%%
feat_data.t = feat_data.t.*1000;
tresh{7} = [tresh{7}.*1000];
%%
[x_tresh] = compute_treshold_xing(feat_data.features,tresh);
%% sum treshold crossing over features
tresh_sig = zeros(1,length(x_tresh{1}));
for it = 1:length(x_tresh)
    tresh_sig = tresh_sig + x_tresh{it};
end;
%% search for speech/vocal units
[tresh_sig] = tresh_sig./tresh{8};
%%
tresh_idx = find(tresh_sig >= 1);
tresh_sig_2tresh = zeros(1,length(tresh_sig));
tresh_sig_2tresh(tresh_idx) =1;

xidx_p= find(diff(tresh_sig_2tresh)==1)+1;
xidx_n = find(diff(tresh_sig_2tresh)==-1)+1;
%% consistency check
if isempty(xidx_p) && ~isempty(xidx_n)
    xidx_p = 1;
end;
if ~isempty(xidx_p) && isempty(xidx_n)
    xidx_n = length(tresh_sig);
end;
if min(xidx_n) < min(xidx_p)
    xidx_n(xidx_n == min(xidx_n)) = [];
end;
if max(xidx_p) > max(xidx_n)
    xidx_n(end+1) = length(feat_data.t);
end;
if isempty(xidx_n)
    xidx_n = length(feat_data.t);
end;
%%
% xidx_n(feat_data.t(xidx_p)<250) = [];
% xidx_p(feat_data.t(xidx_p)<250) = [];
%% syllable segmentation
syl_sig = sign(feat_data.syl_seg.*(abs(feat_data.syl_seg)>=0.25));
syl_idx = find(abs(syl_sig) >=1);
syl_sig_2tresh = zeros(1,length(feat_data.syl_seg));
syl_sig_2tresh(syl_idx) =1;

xidx_p2= find(diff(syl_sig_2tresh)==1)+1;
xidx_n2 = find(diff(syl_sig_2tresh)==-1)+1;
%% consistency check
if min(xidx_n2) < min(xidx_p2)
    xidx_n2(xidx_n2 == min(xidx_n2)) = [];
end;
if max(xidx_p2) > max(xidx_n2)
    xidx_n2(end+1) = max(xidx_p2);
end;
%% consistency check
if length(xidx_p2) == length(xidx_n2)
    for it = 1:length(xidx_p2)
        if xidx_p2(it) >1
            xidx_p2(it)=  xidx_p2(it)-1;
        end;
        if xidx_n2(it) > 1
            xidx_n2(it) = xidx_n2(it)-1;
        end;
    end;
else
    xidx_p2 = [];
    xidx_n2 = [];
end;
%% consistency check
if length(xidx_p) == length(xidx_n)
    %do nothing
else
    xidx_p = [];
    xidx_n = [];
end;
%% consistency check
if (length(xidx_p) == length(xidx_n)) == 0
    error('wrong number of data segments');
end;
%% consistency check
if (length(xidx_p2) == length(xidx_n2)) == 0
    error('wrong number of data segments');
end;
%%
feat_data.finf.xidx_p = NaN;
feat_data.finf.xidx_n = NaN;
feat_data.finf.xidx_p2 = NaN;
feat_data.finf.xidx_n2 = NaN;
feat_data.finf.x_tresh = NaN;
feat_data.finf.features = NaN;
feat_data.finf.tresh = NaN;
feat_data.finf.inf_f = NaN;
feat_data.finf.on_idx = NaN;
feat_data.finf.off_idx = NaN;
feat_data.finf.on_t = NaN;
feat_data.finf.off_t = NaN;
%%
if ~isempty(xidx_p)
    
    %%
    [nsyl] = syllabel_counter(xidx_p,xidx_n,feat_data.t,tresh{7});
    xidx_p = xidx_p(nsyl);
    xidx_n = xidx_n(nsyl);
    %%
    [nsyl2] = syllabel_counter(xidx_p2,xidx_n2,feat_data.t,tresh{7});
    xidx_p2 = xidx_p2(nsyl2);
    xidx_n2 = xidx_n2(nsyl2);
    %%
    if ~isempty(xidx_p)
        
        %%
        if min(xidx_p2) < min(xidx_p)
            xidx_p2(xidx_p2 == min(xidx_p2)) = min(xidx_p);
        end;
        
        if max(xidx_n2) > min(xidx_n)
            xidx_n2(xidx_n2 == max(xidx_n2)) = max(xidx_n);
        end;
        %%
        inf_f = zeros(1,length(x_tresh));
        for it = 1:length(x_tresh)
            if ~isempty(xidx_p2) && ~isempty(xidx_n2)
                inf_f(it) = sum(x_tresh{it}(xidx_p2(1):xidx_n2(1)));
            end;
        end;

        [~,s_idx] = sort(inf_f);
        inf_f = flipud(s_idx')';
        %%
        feat_data.finf.xidx_p = xidx_p;
        feat_data.finf.xidx_n = xidx_n;
        feat_data.finf.xidx_p2 = xidx_p2;
        feat_data.finf.xidx_n2 = xidx_n2;
        feat_data.finf.x_tresh = x_tresh;
        feat_data.finf.features = feat_data.features;
        feat_data.finf.tresh = tresh;
        feat_data.finf.inf_f = inf_f;
        %%
        dt = diff(x_tresh{1});
        nseg = length(find(dt ==1));
        
        ix1 = find(dt==1)+1;
        if isempty(ix1) && x_tresh{1}(1) ==1
            ix1 = 1;
        end;
        ix2 = find(dt==-1)+1;
        
        if length(ix2)>=length(ix1)
            dum = ix2;
            for mt = 1:length(ix2)
                if sign(ix2(mt)-ix1)==-1;
                    dum(mt) = [];
                end;
            end;
            ix2 = dum;
        end;
        
        if length(ix1)>length(ix2)
            ix2(end+1) = length(x_tresh{1});
        end;
        
        seg.idx = cell(length(ix1),1);
        for mt = 1:nseg
            
            seg.idx{mt} = ix1(mt):ix2(mt);
            
        end;
        
        sel_idx = zeros(length(seg.idx));
        k = 0;
        ref = min(xidx_p):max(xidx_n);
        for mt = 1:length(seg.idx)
            if any(ismember(seg.idx{mt},ref))
                k = k+1;
                sel_idx(k) = mt;
            end;
        end;
        sel_idx(k+1:end) = [];
        
        on_idx = min([seg.idx{sel_idx}]);
        off_idx = max([seg.idx{sel_idx}]);
        %         on_idx = min(xidx_p);
        %         off_idx = max(xidx_n);                
        %%
        feat_data.finf.on_idx = on_idx;
        feat_data.finf.off_idx = off_idx;
    end;
end;
%%
if isnan(feat_data.finf.on_idx)
    on = NaN;
    off = NaN;
elseif isempty(feat_data.finf.on_idx)
    on = NaN;
    off = NaN;
else
    on = feat_data.t(feat_data.finf.on_idx);
    off = feat_data.t(feat_data.finf.off_idx);
    %% apply correction factor for length of sliding window
    on = on-10;
    off = off+10;
end;