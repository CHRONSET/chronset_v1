%%
restoredefaultpath;
addpath(genpath('~froux/froux/chronset/'));
%%
%compute_features_SayWhen;
%%
fprintf('features terminated moving on to optimization\n');
%%
mode = 1;
noptim = 250;
nworker = 128;
nparts = 128;
lf = round(nparts/nworker*10)/10;
nworker*lf
%%
%calibration_thresholds_BCN_data_16_June_2016(mode,noptim,nworker,lf);
calibration_thresholds_SayWhen_data_16_June_2016(mode,noptim,nworker,lf,'new');