%% set the MATLAB path
restoredefaultpath;
addpath('C:/fieldtrip-20170619');
%addpath('/home/froux/froux/fieldtrip-20130721/');
ft_defaults;
%%
%path2files = ['/bcbl/data/BEHAVIOUR/BBR_BEH/DATA/Bilingual_Picture_Naming/01_itzal/'];
path2files = ['../demo2/'];
files = dir([path2files,'*.wav']);
%%
% run through the files
for it = 1:1:length(files)
    [wav,Fs] = wavread([path2files,files(it).name]);
    
    wav = wav(:,1)';
    
    dummy = struct;
    dummy.fsample = Fs;
    dummy.trial{1} = wav;
    dummy.time{1} = 0:1/Fs:(length(wav)-1)/Fs;
    dummy.label(1) = {'chan1'};
    
    cfg = [];
    cfg.demean = 'yes';
    cfg.detrend = 'yes';
    cfg.derivative = 'yes';
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 2.5e3;
    
    dummy = ft_preprocessing(cfg,dummy);
    
    cfg = [];
    cfg.detrend = 'no';
    cfg.resamplefs = 5e3;
    
    dummy = ft_resampledata(cfg,dummy);
    
    cfg = [];
    cfg.output = 'pow';
    cfg.method = 'mtmconvol';
    cfg.taper = 'dpss';
    cfg.pad = 'maxperlen';
    cfg.foi = 70:1:dummy.fsample/2;
    cfg.t_ftimwin = 0.075*ones(1,length(cfg.foi));
    cfg.toi = dummy.time{1}(1):0.005:dummy.time{1}(end);
    cfg.tapsmofrq =20;
    
    freq = ft_freqanalysis(cfg,dummy);
    
    dummy2 = struct;
    dummy2.time = freq.freq;
    dummy2.avg = 10*log10(squeeze(freq.powspctrm))';
    dummy2.label = cell(1,length(freq.time));
    for jt = 1:length(dummy2.label)
        dummy2.label(jt) = {['chan',num2str(jt)]};
        dummy2.dimord = 'chan_time';
    end;
    
    cfg = [];
    cfg.output = 'pow';
    cfg.method = 'mtmfft';
    cfg.tapsmofrq = .01;
    cfg.channel = {'all'};
    c = ft_freqanalysis(cfg,dummy2);
    
    %% amplitude
    sf = 10;
    amp = (squeeze(sum(log(freq.powspctrm),2)));
    idx = find(isnan(amp));
    amp(idx) = 0;
    %amp = (amp - mean(amp))./std(amp);
    amp(idx) = NaN;
    amp = conv(gausswin(sf),amp);
    amp = amp(5:end-5);
    
    rge =max(amp)- ((max(amp)-min(amp))/2);
    amp_tresh = amp.*(amp>=rge);
    amp_tresh(isnan(amp_tresh)) = 0;
    amp_tresh(amp_tresh >0) = 1;
    
    ed1 = diff(amp_tresh);
    on_idx1 = find(sign(ed1) ==1);
    off_idx1 = find(sign(ed1) ==-1);
    
    %% harmonic pitch
    sf = 40;
    cep = c.powspctrm;
    hp = max(cep(:,4:end),[],2);
    hp = conv(gausswin(sf),hp);
    
    hp(isnan(hp)) = 0;
    hp = (hp-mean(hp))./std(hp);
    hp = hp(sf/2:end-sf/2);
    hp_tresh = hp.*(hp>0.5);
    hp_tresh(hp_tresh>0)=1;
    
    ed2 = diff(hp_tresh);
    on_idx2 = find(sign(ed2) ==1);
    off_idx2 = find(sign(ed2) ==-1);
    
    %% entropy
    sf = 40;
    df_freq = freq;
    df_freq.powspctrm = diff(df_freq.powspctrm,1,2);
    df_freq.freq = df_freq.freq(2:end);
    
    N = length(freq.freq);
    a = exp(1/N.*(sum(diff(log(freq.powspctrm),1,2),2)));
    b = 1/N.*(sum(freq.powspctrm,2));
    
    a = squeeze(a);
    b = squeeze(b);
    ent = log((a)./(b));
    
    %ent = ent-max(ent);
    
    ent = conv(gausswin(sf),ent);
    ent = ent(sf/2:end-sf/2);
    
    idx1 = find(isnan(ent));
    idx2 = setdiff(1:length(ent),idx1);
    
    z = (ent(idx2) - mean(ent(idx2)))./std(ent(idx2));
    ent(idx2) = z;
    ent(idx1) = 0;
    
    ent_tresh = ent.*(ent<0);
    ent_tresh(ent_tresh>0)=-1;
    
    ed3 = diff(ent_tresh);
    on_idx3 = find(sign(ed3) ==-1);
    off_idx3 = find(sign(ed3) ==1);
    %% get on and off-sets
    idx1 = min(on_idx1):max(off_idx2);
    idx2 = min(on_idx2):max(off_idx2);
    idx3 = min(on_idx3):max(off_idx3);
    
    on_idx = min(intersect(intersect(idx1,idx2),intersect(idx1,idx3)));
    on_t = freq.time(on_idx);
    off_idx = max(intersect(intersect(idx1,idx2),intersect(idx1,idx3)));
    off_t = freq.time(off_idx);
    
    sound_f = wav;
    sound_t = 0:1/Fs:(length(wav)-1)/Fs;
    
    
    %% visualize
    figure;
    subplot(6,1,1:4);
    imagesc(df_freq.time,df_freq.freq(2:end),(squeeze(df_freq.powspctrm)));
    axis xy;colormap gray;caxis([-1e-9 1e-9]);
    ylabel('Frequency [kHz]');
    
    axes('position',[.13 .455 .775 .2]);
    ent2 = ent;
    ent2(ent_tresh>-1) = NaN;
    plot(freq.time,ent2,'m','LineWidth',3);
    axis tight;axis off;
    
    axes('position',[.13 .455 .775 .2]);
    amp2 = amp;
    amp2(amp_tresh<1) = NaN;
    plot(freq.time,amp2,'y','LineWidth',3);
    axis tight;axis off;
    
    axes('position',[.13 .575 .775 .2]);
    hp2 = hp;
    hp2(hp_tresh<1) = NaN;
    plot(freq.time,hp2,'g','LineWidth',3);
    axis tight;axis off;
    
    subplot(6,1,5:6);
    hold on;
    plot(0:1/Fs:(length(wav)-1)/Fs,wav);
    for jt = 1:length(on_idx)
        plot(freq.time(on_idx(jt)),0,'g*');
        plot(freq.time(off_idx(jt)),0,'g*');
    end;
    plot([freq.time(min(on_idx)) freq.time(min(on_idx))],[min(wav) max(wav)],'r');
    plot([freq.time(max(off_idx)) freq.time(max(off_idx))],[min(wav) max(wav)],'r');
    plot([freq.time(min(on_idx)) freq.time(max(off_idx))],[min(wav) max(wav)],'r');
    plot([freq.time(min(on_idx)) freq.time(max(off_idx))],[max(wav) min(wav)],'r');
    axis tight;
    xlabel('Time [s]');
    ylabel('Amplitude [a.u.]');
    
        %%
        flag = 0;
        while flag <1
            [x,y] = ginput(2);
            time = 0:1/Fs:(length(wav)-1)/Fs;
            idx = find(time >= x(1) & time <= x(2));
            
            sel_sig = wav(idx);
            
            soundsc(sel_sig,Fs);
            
            [s] = input('do you want to keep parameters? [y/n]','s');
            if strcmp(s,'y')
                flag = 1;
            end;
        end;
    
end;
