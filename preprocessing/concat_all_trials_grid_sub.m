function data_all_trials = concat_all_trials_grid_sub(Trials,experiment,stop_condition,...
         driveNameAnalyze,MONKEYDIR)
% this function pool all recordings from all trials for each day and return
% in the from electrode * time
% by all data for each trial, we mean from StartOn of the first trial to
% ReachStop of the last trial
%Inputs
%   Trials:the regular trial structure
%   stop_condition: as per definition of a trial for either "ReachStop" or
%   "TargOn"
%   driveNameAnalyze: 'LM1_ECOG_3' or 'LM1_SC32'
%Outputs
%   data_all_trials: concatenated data of the form [electrode time]

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
trial_start_ReachEnd = zeros(length(Trials),1);

data_all_trials = [];
for ti = 1:length(Trials)
    trial_temp = Trials(ti);
  
    trialInfo.trig = 'StartOn';
    trialInfo.tBefore = 0;
    
    %check the stop condition
    if strcmp(stop_condition, 'TargsOn')
     trialInfo.tAfter  = trial_temp.TargsOn;
    elseif strcmp(stop_condition, 'ReachStop')
        trialInfo.tAfter  = trial_temp.ReachStop;
    else
        error('Unrecognized Stop Condition: ' + stop_condition)
    end
    
    trial_start_ReachEnd(ti) = trialInfo.tAfter;
    
    data = trialLfp(trial_temp,driveNameAnalyze,dataParam.eList,[],...
           trialInfo.trig,[trialInfo.tBefore trialInfo.tAfter],trialInfo.lfpType,MONKEYDIR);
    
    %do grid subtraction
    data = ref_grid_sub(data);
       
    %concatenate all trials
    data_all_trials = [data_all_trials data];
end


% data = trialLfp(Trials,driveNameAnalyze,dataParam.eList,[],trialInfo.trig,[],trialInfo.lfpType,MONKEYDIR);

% %second way
% trial_temp = Trials(1); %from first StartOn to last ReachStop
% trialInfo.trig = 'StartOn';
% trialInfo.tBefore = 0;
% trialInfo.tAfter = Trials(end).StartOn + Trials(end).ReachStop;
% 
% data_all_trials = trialLfp(trial_temp,driveNameAnalyze,dataParam.eList,[],...
%             trialInfo.trig,[trialInfo.tBefore trialInfo.tAfter],trialInfo.lfpType,MONKEYDIR);

disp('Finshed concatenating all trials')
end


%code archive, this is to figure out 
