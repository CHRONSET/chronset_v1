function [on] = chronset_textUI()

% chronset_textUI(params)
%[Inputs]: params is a strcuture with different fields
%
%	mandatory inputs:
%	-----------------
%	params.input_file_name 	= filename of the audiofile to be processed
%
%
%	optional inputs:
%	----------------
%	params.savepath 	= path to where the output txt file will be saved
%
%
%[Output]: this function returns the onset time of speech in 2 formats:
%          1) the onset is printed out on the screen
%          2) a txt file is saved to disc to either the default path or a
%          path specified by the user
%              
%
% F.Roux, University of Birmingham, Jul 2016
	

%% set the path to the default chronset folder
restoredefaultpath;
if isunix
	chrondir = dir('/usr/local/chronset/');
	if isempty(chrondir)
		error('search fro chronset in /usr/local/ failed. please consult README.txt');
	end;
    chrondir = '/usr/local/chronset/';
end;

if ispc
    chrondir = dir('C:\Program Files\');
    if isempty(chrondir)
        error('search fro chronset in /usr/local/ failed. please consult README.txt');
    end;
    chrondir = 'C:\Program Files\';
end;

addpath(genpath(chrondir));

%% get the path and filename of the input wav-file

[FILENAME, LOAD_PATHNAME, ~] = uigetfile(chrondir,'Welcome to Chronset. pick the files you wish to load','MultiSelect', 'on');
[SAVE_PATHNAME] = uigetdir(chrondir,'Thank you! Now pick a directory to save your output');

if iscell(FILENAME) <2
    dum{1} = FILENAME;
    FILENAME = dum;
    clear dum;
end;

%% load the precomputed threshods

load([chrondir,filesep,'thresholds',filesep,'greedy_optim_thresholds_BCN_final.mat']);

[thresh] = chronset_extract_thresholds(optim_data);

%% reads in the wavfile
tt = tic;
for it = 1:length(FILENAME)
    in = struct;
    fprintf('computing file\n');
    [in.wav,in.FS] = wavread2([LOAD_PATHNAME,FILENAME{it}]);
    
    %% compute speech features
    
    [feat_data] = compute_feat_data([],in);
    
    %% detect speech onset
    [on(it)] = detect_speech_on_and_offset(feat_data,[thresh' {0.035} {4} {0.25}]);
    
    %% write output textfile
    fid = fopen([SAVE_PATHNAME,filesep,FILENAME{it},'_onset.txt'],'w+');
    fprintf(fid,'%s',FILENAME{it});
    fprintf(fid,'\t');
    fprintf(fid,'%s',num2str(on));
    fprintf(fid,'\n\r');
    
end;
tt = toc(tt);

fprintf('Speech onset detection complete\n');
fprintf(['Total time required:',num2str(round(tt*1e2)/1e2),' seconds \n']);
fprintf('Thank you for using CHROSET!!!\n');
fprintf('Please leave us your feedback at infox@chron.net\n');
