chronset_batch('./data/wav_file_example/','./data/results.txt')


chronset_batch('/bcbl/home/home_a-f/barmstrong/agnesa/agnesa_nans/','/bcbl/home/home_a-f/barmstrong/agnesa/Noise20.txt')


chronset_batch('/bcbl/home/home_a-f/barmstrong/agnesa/UCDavis/','/bcbl/home/home_a-f/barmstrong/agnesa/davis.txt')




%%%
it = 1;
LOAD_PATHNAME = '/bcbl/home/home_a-f/barmstrong/agnesa/UCDavis/';
FILENAME{it} = 'ISPC_B_V2-23-3Slide1-105.wav';
load([chrondir,filesep,'thresholds',filesep,'greedy_optim_thresholds_BCN_final.mat']);

[thresh] = chronset_extract_thresholds(optim_data);


in = struct;
[in.wav,in.FS] = wavread2([LOAD_PATHNAME,FILENAME{it}]);
in.wav = in.wav(:,1);
%replace completely empty parts of a recording with low amount of noise

if fillEmpty == 1
    %lock random noise values for replication purposes.
    rng(1);
    
    sig = in.wav;
    
    if ~isempty(sig(sig~=0))
        A = median(sig(sig~=0));
    else
        A = 0.01;
    end
    
    if ~isempty(sig(sig~=0))
        %20 percentile seems enough noise to avoid NaNs due to
        %singularity
        qt = quantile(sig(sig~=0),0.2);
    else
        qt = 0.01;
    end
    
    zlx = find(sig==0);
    
    sig(zlx) = A +qt*randn(1,length(zlx));
    
    in.wav = sig;
end


%% compute speech features

[feat_data] = compute_feat_data([],in);

feat_data = feat_data;
tresh = [thresh' {0.035} {4} {0.25}];

%[on(it)] = detect_speech_on_and_offset(feat_data,[thresh' {0.035} {4} {0.25}]);