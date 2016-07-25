function [features,syl_seg] = compute_Speech_features(feat_data)
%% indexes for frequencies of interest
freq_sel_idx = 1:length(feat_data.f);%find(feat_data.f > 50);
%% amplitude in dB - feature #1
%[amp] = 20*log10(feat_data.s+1);
[amp] = 20*log10(feat_data.s+1+eps);
[amp] = sum(amp,2);
%% wiener entropy - feature # 2
N = length(feat_data.f);
a = exp(1/N*sum(log(feat_data.s(:,freq_sel_idx)),2));
b = 1/N*sum(feat_data.s,2);
[we] = log(a)-log(b);
%% change over time and frequency - feature #3
x = (feat_data.ds(1,:,:).^2.*feat_data.ds(2,:,:));
x = max(x,[],3);
T = max(x)/100*5;
[dfdt] = x.*(x>T); 
%% amplitude modulation (AM) - feature # 4
[am] = squeeze(sum(feat_data.ds(1,:,freq_sel_idx).^2,3));%overall time derivative spower -> positive at beginning of sounds & negative at the end of sounds
%% frequency modulation (FM) - feature # 5
fm = atan(max(feat_data.ds(1,:,:).^2,[],3)./max(feat_data.ds(2,:,:).^2,[],3));
%% goodness of pitch - feature #6
C = log(abs(fft(squeeze(feat_data.ds(2,:,:)),[],2)).^2)+log(abs(fft(squeeze(feat_data.ds(1,:,:)),[],2)).^2);
[goP] = max(C,[],2);
%% prepare feature data for output
features =cell(6,1);
features{1} = amp;%amplitude
features{2} = we;%entropy
features{3} = dfdt;
features{4} = am;%amplitude modulation
features{5} = fm;%frequency modulation
features{6} = goP;%Cepstrum

%check consistency of dims
for it = 1:length(features)
    
    if size(features{it},2) > size(features{it},1)
        features{it} = features{it}';
    end;
    
end;
%% syllable segmentation
% z-score for syllable segmentation
x = conv(diff(features{4}),hanning(54));
x = x(28:end-27);

x = x-mean(x);
x = x./std(x);
syl_seg = x;
%% apply transformation to AM and FM
features{4} = diff(cumsum(abs(features{4}))./sum(abs(features{4}))); % 
features{4} = [features{4}(1);features{4}];
