%first we pull trials from different day

%develope the code pipeline to process lfp data
global MONKEYDIR
MONKEYDIR = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\data';
addpath(genpath([MONKEYDIR '/m/']))

% look at the sessions
drive_base = 'LM1_ECOG';
driveSessions = makeDriveDatabase(drive_base,{'180328','180415'});
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

trialInfo.tBefore = -1.5e3;
trialInfo.tAfter  = 1e3;
trialInfo.timeVec = linspace(trialInfo.tBefore,trialInfo.tAfter, trialInfo.tAfter - trialInfo.tBefore);
trialInfo.timeReachStart = abs(trialInfo.tBefore);
trialInfo.trig = 'ReachStart';
trialInfo.lfpType = 'lfp';
trialInfo.Fs_lfp = 1000;
trialInfo.badECoGs = [47 59 163];
trialInfo.goodECoGs = setdiff(1:211,trialInfo.badECoGs);

%define neural feature
featureROI  = [-50 200]; %in ms

%we are going to look at all electrode,  211 ECOG and 32 SC32 electrodes
trialInfo.proc_electrodes = 1:243;

%use constant target defination
selectedfile = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
load(selectedfile)

%we load the regular trials
SessAnalyze = Sessions(useSess,:);
nSess = size(SessAnalyze,1);

%%


for iD=1:nSess
    %load the trials and get the depth information
    trFN = [MONKEYDIR '/' DayAnalyze{iD} '/mat/Trials.mat'];
    load(trFN,'Trials')
    trialInfo.sessName =  DayAnalyze{iD};
    trialInfo.depthProfile = (Trials(1).Depth{1,2})'; % in micron
    depthProfile_all(iD,:,:) = trialInfo.depthProfile;
    
    %filter out trials
    [trialInfo,Trials] = filter_trials(trialInfo,Trials);
    
    %first load an experiment trial
    expFN = [MONKEYDIR '/' DayAnalyze{iD} '/' Trials(1).Rec '/rec' Trials(1).Rec '.experiment.mat'];
    load(expFN,'experiment')
    
    %load data
    [trLfpData,trialInfo] = load_data(Trials, experiment,trialInfo,MONKEYDIR);
    
    %set these electrodes to nan
    trLfpData(:,trialInfo.badECoGs,:) = nan;
    
    %concatenate the trials
    trLfpData_cat = compress_dim(trLfpData,2); %2 preserve electrode dim
    
    %calculate the low frequency components
    r_electrodes = corrcoef(trLfpData_cat);
    
    r_electrodes_all(iD,:,:) = r_electrodes;
    
    disp(['Finished ' DayAnalyze(iD)])
        
    
end

%%
ECOG_offset = 211;

electrodes = [6 7 17 31];

for i = 1:length(electrodes)
    ei = electrodes(i);
    r_ei = squeeze(nanmean(r_electrodes_all(:,1:ECOG_offset,ei + ECOG_offset),2));
    depth_Profile_ei =depthProfile_all(:,ei);
    [B,I] = sort(depth_Profile_ei);
    I(B<=0) =[];
    B(B<=0) = [];
    scatter(B,r_ei(I));
    hold on
end

xlabel('Lowered Distance (um)')
ylabel('Coefficient')
title('Mean correlation coefficient')
hold off
legend(string(electrodes ))
