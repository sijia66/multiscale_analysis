%% load some example recordings
global MONKEYDIR
%MONKEYDIR = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\data';
MONKEYDIR = 'C:\Users\Si Jia\OneDrive - UW\projects\Brain EEG\data'
addpath(genpath([MONKEYDIR '/m/']))

% look at the sessions
drive_base = 'LM1_ECOG';
driveSessions = makeDriveDatabase(drive_base,{'180328','180328'});
Sessions = cat(1, driveSessions{:});
driveNames = Sessions(:,3);
driveNameAnalyze = {'LM1_ECOG_3'};
useSess          = ismember(driveNames, driveNameAnalyze);
DayAnalyze = Sessions(useSess,1);

% define processing information
trialInfo.sequence_random = 1;

%yes, we are filtering out the long acq time
trialInfo.filterAcq = 1;
trialInfo.filterAcq_time = 5000; %ms

trialInfo.tBefore = -1e3; %ms
trialInfo.tAfter  = 1e3;
trialInfo.timeVec = linspace(trialInfo.tBefore,trialInfo.tAfter, trialInfo.tAfter - trialInfo.tBefore);
trialInfo.timeReachStart = abs(trialInfo.tBefore);
trialInfo.trig = 'ReachStart';
trialInfo.lfpType = 'lfp';
trialInfo.Fs_lfp = 1000;
trialInfo.badECoGs = [47 59 163];
trialInfo.ECoG_offset = 211;
trialInfo.goodECoGs = setdiff(1:211,trialInfo.badECoGs);


%we are going to look at all electrode,  211 ECOG and 32 SC32 electrodes
trialInfo.proc_electrodes = 1:243;

%use constant target defination
selectedfile = 'C:\Users\Si Jia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
load(selectedfile)

%we load the regular trials
SessAnalyze = Sessions(useSess,:);
nSess = size(SessAnalyze,1);

for iD=1:nSess
    %load the trials and get the depth information
    trFN = [MONKEYDIR '/' DayAnalyze{iD} '/mat/Trials.mat'];
    load(trFN,'Trials')
    trialInfo.sessName =  DayAnalyze{iD};
    trialInfo.depthProfile = (Trials(1).Depth{1,2})'; % in micron
    trialInfo.goodSC32 = (find( trialInfo.depthProfile > 0))'; 
    trialInfo.goodE = [trialInfo.goodECoGs trialInfo.goodSC32+trialInfo.ECoG_offset];
    trialInfo.badE = setdiff(trialInfo.proc_electrodes,trialInfo.goodE);
    
    %filter out trials
    [trialInfo,Trials] = filter_trials(trialInfo,Trials);
    
    %first load an experiment trial
    expFN = [MONKEYDIR '/' DayAnalyze{iD} '/' Trials(1).Rec '/rec' Trials(1).Rec '.experiment.mat'];
    load(expFN,'experiment')
    
    %load data
    [trLfpData,trialInfo] = load_data(Trials, experiment,trialInfo,MONKEYDIR);
    
    trLfpData(:,trialInfo.badECoGs,:) = nan;
    
        %do target sorting
    [Pos, targids] = calcSeqReachTarg(Trials);
    Pos_Seq_1 = squeeze(Pos(:,2,:));
    trial_TargetAssigned = assignTaskNumber(Pos_Seq_1_unique,Pos_Seq_1);
end
%% build the feature set

S_Z_all_electrodes = [];
%feature_electrode_location = [];
for i = 1:length(trialInfo.proc_electrodes)
    example_electrode = trialInfo.proc_electrodes(i);
    
    dat = squeeze(trLfpData(:,example_electrode,:));
    data = dat'; % times * trials
    
    %calculate the  spectrum
    [S,t,f] = mtspecgramc( data, movingwin, params );
    t = t + trialInfo.tBefore/1000; %shift relative to movement start
    S_Z = zscore(S,[],1);
    
    %save this to a feature vector
    S_Z_all_electrodes(i,:,:,:) = ...
        S_Z;
    %notify the real status of the code
    if mod(example_electrode,10) == 0
        disp(['finished: ',num2str(example_electrode)])
    end
end
%% decode the target directions the smart way
%only select the electrodes that predict the target direction

pred_SC32 = intersect(find(feature_electrode_location(:,1) < 0 )...
                         ,trialInfo.goodSC32 + trialInfo.ECoG_offset)...
            - trialInfo.ECoG_offset;
        
ind_SC32 = intersect(find(feature_electrode_location(:,1) > 0 )...
                         ,trialInfo.goodSC32 + trialInfo.ECoG_offset)...
            - trialInfo.ECoG_offset;
        
pred_ECoG = intersect(find(feature_electrode_location(:,1) < 0 )...
                         ,trialInfo.goodECoGs);
                     
ind_ECoG  = intersect(find(feature_electrode_location(:,1) > 0 )...
                         ,trialInfo.goodECoGs);
                     
                     
%get frequency features for SC32
for i = 1:length(ind_SC32)
    ei = ind_SC32(i) + trialInfo.ECoG_offset;
    time_ind = find(feature_electrode_location(ei,1) == t);
    freq_ind = find(feature_electrode_location(ei,2) == f);
    freq_features(i,:) = squeeze(S_Z_all_electrodes(ei,time_ind,freq_ind,:));
end

[cat_hat] = runLeaveOneOutClassification(freq_features', trial_TargetAssigned)

accuracy = sum(cat_hat == trial_TargetAssigned) / length(trial_TargetAssigned)

