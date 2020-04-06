function data_all_trials = concat_all_electrodes(trial_temp,experiment, driveNameAnalyze,MONKEYDIR)
% this function pool all recordings from one trial for all electrodes and return
% in the from  [time]
% by all data for each trial, we mean from StartOn to ReachStop of the first trial to
%Inputs
%   trial_temp:the regular trial structure
%   driveNameAnalyze: 'LM1_ECOG_3' or 'LM1_SC32'
%Outputs
%   data_all_trials: concatenated data of the form [1 time]

% this needs to be changed

if strcmp(driveNameAnalyze,'LM1_ECOG_3')
    drive_idx = 1;
elseif strcmp(driveNameAnalyze,'LM1_SC32_1')
    drive_idx = 2;
else
    error('Unrecognized Drive Name: ' + driveNameAnalyze)
end


% trialInfo.tBefore = -2e3;
% trialInfo.tAfter  = 1e3;
trialInfo.lfpType = 'lfp';
trialInfo.Fs_lfp = 1000;

%inspect the trials and figure out the shortest period between startOn and
%reachStop
dataParam.eList = [experiment.hardware.microdrive(drive_idx).electrodes.channelid];


data_all_trials = [];
% for ti = 1:length(Trials)
%     trial_temp = Trials(ti);
%   
%     trialInfo.trig = 'StartOn';
%     trialInfo.tBefore = 0;
%     trialInfo.tAfter  = trial_temp.TargsOn;
%     
%     trial_start_ReachEnd(ti) = trialInfo.tAfter;
%     
%     data = trialLfp(trial_temp,driveNameAnalyze,dataParam.eList,[],...
%            trialInfo.trig,[trialInfo.tBefore trialInfo.tAfter],trialInfo.lfpType,MONKEYDIR);
%     %concatenate all trials
%     data_all_trials = [data_all_trials data];
% end


% data = trialLfp(Trials,driveNameAnalyze,dataParam.eList,[],trialInfo.trig,[],trialInfo.lfpType,MONKEYDIR);

%second way
trialInfo.trig = 'StartOn';
trialInfo.tBefore = 0;
trialInfo.tAfter = trial_temp(1).ReachStop;

data_temp = trialLfp(trial_temp,driveNameAnalyze,dataParam.eList,[],...
            trialInfo.trig,[trialInfo.tBefore trialInfo.tAfter],trialInfo.lfpType,MONKEYDIR);

data_all_trials = reshape(data_temp',1,[]);

disp('Finshed concatenating all trials')
end


%code archive, this is to figure out 
