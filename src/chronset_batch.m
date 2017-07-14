%Usage after compile:
% % ./run_chronset_batch.sh <MCRDEPLOY_ROOT> <dir with input wavs> <outfile> 
%
% e.g.,
% ./run_chronset_batch.sh /opt/matlab/MATLAB_Compiler_Runtime/v80/ ~/barmstrong/Chronset_input/ chronOut_test.txt



function chronset_batch(input_folder,output_file)

%write any errors to output_file rather than dying quietly

fillEmpty = 1;

figSaveDir = '/bcbl/home/home_a-f/barmstrong/agnesa/figs/';

plotSpeechFeatures = 1;

nWorkers = 12; % number of parallel workers

dropLast25ms = 0; %only works if fillEmpty is also true


%% Load Optimized Thresholds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['.',filesep,'thresholds',filesep,'greedy_optim_thresholds_BCN_final.mat']);
%%
[thresh] = chronset_extract_thresholds(optim_data);



%% Process Individual Wavs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    fileList = dir(input_folder);
    nf = length(fileList);
    rts = cell(nf,1);

    matlabpool(nWorkers);

    %start at ind = 3 to avoid . and .. on the search path
    %parfor i = 3:nf
    for i = 3:nf
        %Disp is now within individual file call, later.
        %disp(['File name being processed: ' fileList(i).name]);

        in = struct;
        try
            [in.wav,in.FS] = wavread2([input_folder '/' fileList(i).name]);
        catch ME
            [in.wav,in.FS] = audioread([input_folder '/' fileList(i).name]);
        end
    
        if ~isempty(in.wav)
            %[in.wav,in.FS] = wavread2([input_folder '/' fileList(i).name]);
            in.wav = in.wav(:,1);

            %in.wav = in.wav(:,1);
        %replace completely empty parts of a recording with low amount of noise

        if fillEmpty == 1
              sig = in.wav;


            %drop last 25 ms from files if you using the fMRI denoising script
            if dropLast25ms ==1
                %disp(length(sig));
                %disp(in.FS);
                %disp('dropping last 25 ms');
                %sig(end-250:end) = 0;
                %sig(1:250) = 0;

                sig(1:round(0.05*in.FS)) = 0;
                sig(round(end-0.05*in.FS):end) = 0;
                %sig(0.1*in.FS:end-0.1*in.FS) = 0;
            end

            %lock random noise values for replication purposes.  
            rng(1);



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



            try
                [feat_data] = compute_feat_data([],in);
                [on] = detect_speech_on_and_offset(feat_data,[thresh' {0.035} {4} {0.25}]);

                rts(i) = {[fileList(i).name, '	', num2str(round(on))]};


                if plotSpeechFeatures == 1
                    plot_speech_features(feat_data, thresh,figSaveDir,[fileList(i).name,'.pdf']);
                end

            catch ME
                disp('Error processing file.  Chronset will try to continue crushing RTs...');
                disp(getReport(ME));
                rts(i) = {[fileList(i).name, '	', 'ERROR_WAVBAD?']};
            end
            
        else
            rts(i) = {[fileList(i).name, '	', 'EMPTYFILE']};
        end
        
    end;

    matlabpool close;
catch ME
    
    msgString = getReport(ME);
    rts = cell(3,1);
    rts(3) = {msgString};
    
end;

%% Saving Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Writing results to file:');

disp(output_file);


fid = fopen(output_file,'w+');

nrows = size(rts);
formatSpec = '%s\n';

%b/c ind = 3 used in processing individual wavs
for row = 3:nrows
	fprintf(fid,formatSpec,rts{row,:});
end

fclose(fid);
