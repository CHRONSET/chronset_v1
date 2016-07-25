function [features] = compute_Speech_features2(feat_data)
%%
fdt = squeeze(feat_data.ds(1,:,:));
fdf = squeeze(feat_data.ds(2,:,:));
%% amplitude in dB - feature #1
[amp] = 20*log10(feat_data.s+1);
[amp] = max(amp,[],2);
%% wiener entropy - feature # 2
N = length(feat_data.f);
a = exp(1/N*sum(log(feat_data.s),2));
b = 1/N*sum(feat_data.s,2);
[we] = log(a)-log(b);
%% change over time and frequency - feature #3
x = (fdt.^2+fdf.^2);
x = max(x,[],2);
T = max(x)/100*5;
[dfdt] = x.*(x>T); 
%% amplitude modulation (AM) - feature # 4
[am] = sum(fdt.^2,2);%overall time derivative spower -> positive at beginning of sounds & negative at the end of sounds
[am] = abs(am);
%% frequency modulation (FM) - feature # 5
[fm] = atan(max(fdt.^2,[],2)./(max(fdf.^2,[],2)));
%% goodness of pitch - feature #6
C = log(abs(fft(fdf,[],2)).^2);
[goP] = mean(C,2);
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
%% apply transformation to AM and FM
features{4} = diff(cumsum(abs(features{4}))./sum(abs(features{4}))); % 
features{4} = [features{4}(1);features{4}];