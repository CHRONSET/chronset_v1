function [on,off,feat_data] =  detect_speech_on_and_offset(feat_data,tresh,mS)
% Use as: [on,off,feat_data] = detect_speech_on_and_offset(path2file,file_name,tresh)
%
% Function returns the on- and off-sets of vocal signals based on feature data (units are in ms).
% Input: 
% -in_data: structure with fields in_data.wav and in_data.FS as returned by wavread2
% - tresh = 1x6 cell array that contains the treshold values 
% Output:                          
%- on:  onset time,                                  
%- off: offset time                               
%- feat_data: structure array with feature data derived from spectral analysis
%
% default threshold parameters:
% tresh = cell(length(feat_data.features),1);
% tresh{1} = 5e-3;%amplitude
% tresh{2} = .55;%wiener entropy
% tresh{3} = 5e-3;%spectral change
% tresh{4} = 5e-3;% amplitude modulation
% tresh{5} = .6;%frequency modulation
% tresh{6} = .25;%goodness of pitch

% Author: F. Roux 
% Basque Center on Cognition, Brain & Language (BCBL),
% Donostia San-Sebastian, September 2015
%% convert time units from s to ms
feat_data.t = feat_data.t.*1000;
%% binarize feature vectors
[x_tresh] = compute_treshold_xing(feat_data.features,tresh);
%% compute significance statistic of thresholded feature data 
tresh_sig = zeros(1,length(x_tresh{1}));
for it = 1:length(x_tresh)%
    tresh_sig = tresh_sig + x_tresh{it};
end;
%% binarize significance statistic
[tresh_sig] = tresh_sig./4;% at least 3 features must cross threshold
tresh_idx = find(tresh_sig >= 1);%search for significant time stamps

tresh_sig_2tresh = zeros(1,length(tresh_sig));
tresh_sig_2tresh(tresh_idx) =1;% set significant time stamps to 1

xidx_p= find(diff(tresh_sig_2tresh)==1)+1;% indexes of speech beginings
xidx_n = find(diff(tresh_sig_2tresh)==-1)+1;% indexes of speech endings
%% consistency check
if min(xidx_n) < min(xidx_p)
    xidx_n(xidx_n == min(xidx_n)) = [];%xidx_p = [1 xidx_p];%this makes sure no speech gets carried over
end;
if max(xidx_p) > max(xidx_n)
    xidx_n(end+1) = length(feat_data.t);%just for consistency in case speech extends outside of recording
end;
if isempty(xidx_n) && ~isempty(xidx_p)
    xidx_n = length(feat_data.t);
end;
%% flag to filter freak onsets
% xidx_n(feat_data.t(xidx_p)<.25) = [];
% xidx_p(feat_data.t(xidx_p)<.25) = [];
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
%% intialize data for saving results
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
%% do onset detection (aka crush some latencies)
if ~isempty(xidx_p)
    
    %% check sylable unit duration (must be above 35 ms)
    [nsyl] = syllabel_counter(xidx_p,xidx_n,feat_data.t);
    xidx_p = xidx_p(nsyl);
    xidx_n = xidx_n(nsyl);
    %%
    if ~isempty(xidx_p)
        %% rank features according to vote contribution
        inf_f = NaN(1,length(x_tresh));
        c = 0;
        sel_idx = min(xidx_p):max(xidx_n);
        if ~isempty(mS)
            for it = 1:length(x_tresh)%
                c = c+1;%ranking feature distance from manual score
                ix = min(find(x_tresh{it} == 1));
                if ix ~=1
                    ix = sel_idx(find(x_tresh{it}(sel_idx)==1));
                    if ~isempty(ix)
                        inf_f(c) = abs(feat_data.t(ix(1))-mS);
                    end;
                end;
            end;
            
            v = min(find(inf_f == min(inf_f)));
        else
            for it = 1:length(x_tresh)%
                c = c+1;%ranking feature distance from total sign samples
                inf_f(c) = sum(x_tresh{it}(min(xidx_p):max(xidx_n)));
            end;
            
            [~,s_idx] = sort(inf_f);%do the ranking
            %v = s_idx(1);
            v = 1;
        end;
        %% save info
        feat_data.finf.xidx_p = xidx_p;
        feat_data.finf.xidx_n = xidx_n;
        feat_data.finf.x_tresh = x_tresh;
        feat_data.finf.features = feat_data.features;
        feat_data.finf.tresh = tresh;
        feat_data.finf.inf_f = inf_f;
        
        %         %
        %         % readout those segments where amplitude and at least 2 other non-amplitude features were simultaneously above threshold
        %         %if inf_f(1) ==1
%         dt = diff(x_tresh{v});
%         nseg = length(find(dt ==1));
%         
%         ix1 = find(dt==1)+1;
%         if isempty(ix1) && x_tresh{v}(1) ==1
%             ix1 = 1;
%         end;
%         ix2 = find(dt==-1)+1;
%         
%         if length(ix2)>=length(ix1)
%             dum = ix2;
%             for mt = 1:length(ix2)
%                 if sign(ix2(mt)-ix1)==-1;
%                     dum(mt) = [];
%                 end;
%             end;
%             ix2 = dum;
%         end;
%         
%         if length(ix1)>length(ix2)
%             ix2(end+1) = length(x_tresh{v});
%         end;
%         
%         seg.idx = cell(length(ix1),1);
%         for mt = 1:nseg
%             
%             seg.idx{mt} = ix1(mt):ix2(mt);
%             
%         end;
%         
%         %compute onset based on amplitude feature
%         sel_idx = zeros(length(seg.idx),1);
%         k = 0;
%         ref = min(xidx_p):max(xidx_n);
%         for mt = 1:length(seg.idx)
%             if any(ismember(seg.idx{mt},ref))
%                 k = k+1;
%                 sel_idx(k) = mt;
%             end;
%         end;
%         sel_idx(k+1:end) = [];
%         
%         %get onset
%         on_idx = min([seg.idx{sel_idx}]);
%         %get offset
%         off_idx = max([seg.idx{sel_idx}]);
        %else
        %% get on and offset indexes
        ix = find(x_tresh{v}==1);
        on_idx = ix(1);%min(xidx_p);
        off_idx = ix(end);%max(xidx_n);
        %end;
        %         % apply correction factor for half length of sliding window
        %         cf = floor((0.01)/2*(1/((feat_data.t(2)-feat_data.t(1))./1000)));
        %         on_idx = on_idx+cf;%add cf to onset index
        %         off_idx = off_idx-cf;%remove cf from offset index
        %% save onset index
        feat_data.finf.on_idx = on_idx;
        feat_data.finf.off_idx = off_idx;
    end;
end;
%%
if isnan(feat_data.finf.xidx_p) 
    on = NaN;
    off = NaN;
elseif isempty(feat_data.finf.on_idx)
    on = NaN;
    off = NaN;
else
    on = feat_data.t(feat_data.finf.on_idx);
    off = feat_data.t(feat_data.finf.off_idx);
end;